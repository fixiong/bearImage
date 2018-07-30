#include "possion_stiching.h"
#include "possion_solver.h"

#include <assert.h>
#include <cmath>

#include <algorithm>
#include <vector>

#include "utility.hpp"

using namespace bear;


static void x_d(const PImage &img)
{
	for (int y = 0; y < height(img); ++y)
	{
		unsigned short * row = (unsigned short *)scanline(img, y);
		for (int x = width(img) - 1; x > 0; --x)
		{
			row[x] = row[x] - row[x - 1] + 32768;
		}

		row[0] = 32768;
	}
}


static void y_d(const PImage &img)
{
	for (int y = height(img) - 1; y >= 0; --y)
	{
		unsigned short * row_1 = (unsigned short *)scanline_bound(img, y - 1);
		unsigned short * row = (unsigned short *)scanline(img, y);
		for (int x = 0; x < width(img); ++x)
		{
			row[x] = row[x] - row_1[x] + 32768;
		}
	}
}


template<typename T>
struct _ZERO_PS {};

template<>
struct _ZERO_PS<float>
{
	constexpr static float run()
	{
		return 0.0f;
	}

	static float to_unit(int v)
	{
		return (float)v;
	}

	static void to_char(unsigned char &dst,float v)
	{
		dst = limiteU8((int)(dst - v + 0.5f));
	}
};

template<>
struct _ZERO_PS<unsigned short>
{
	constexpr static unsigned short run()
	{
		return 32768;
	}

	static unsigned short to_unit(int v)
	{
		return (v << 7) + 32768;
	}

	static void to_char(unsigned char &dst, unsigned short v)
	{
		dst = limiteU8(dst - round_shift<7>((int)v - 32768));
	}
};


template<typename Unit>
static void x_d(
	const PImage &dx,
	const PImage &src,
	unsigned int ch,
	const PImage &mask,
	int max_grandient)
{
	typedef _ZERO_PS<Unit> ZP;

	int mg = max_grandient & 255;
	for (int y = 0; y < height(src); ++y)
	{
		Unit * drow = (Unit *)scanline(dx, y);
		unsigned char * srow = (unsigned char *)scanline(src, y);
		unsigned char * mrow = (unsigned char *)scanline(mask, y);

		int bpp = src.n_channel;

		drow[0] = ZP::run();
		srow += bpp;

		for (int x = 1; x < width(src); ++x)
		{
			bool f1 = 0 != mrow[x - 1];
			bool f2 = 0 != mrow[x];

			if (f1 == f2)
			{
				drow[x] = ZP::run();
			}
			else
			{
				int g = srow[ch] - (srow - bpp)[ch];

				g = std::min(g, mg);
				g = std::max(g, -mg);

				drow[x] = ZP::to_unit(g);
			}

			srow += bpp;
		}
	}
}

template<typename Unit>
static void y_d(
	const PImage &dy,
	const PImage &src,
	unsigned int ch,
	const PImage &mask,
	int max_grandient)
{
	typedef _ZERO_PS<Unit> ZP;

	int mg = max_grandient & 255;

	Unit * drow = (Unit *)scanline(dy, 0);
	for (int x = 0; x < width(src); ++x)
	{
		drow[x] = ZP::run();
	}

	for (int y = 1; y < height(src); ++y)
	{
		Unit * drow = (Unit *)scanline(dy, y);
		unsigned char * srow = (unsigned char *)scanline(src, y);
		unsigned char * mrow = (unsigned char *)scanline(mask, y);

		unsigned char * srow_ = (unsigned char *)scanline(src, y - 1);
		unsigned char * mrow_ = (unsigned char *)scanline(mask, y - 1);

		int bpp = src.n_channel;

		for (int x = 0; x < width(src); ++x)
		{
			bool f1 = 0 != mrow_[x];
			bool f2 = 0 != mrow[x];

			if (f1 == f2)
			{
				drow[x] = ZP::run();
			}
			else
			{
				int g = srow[ch] - srow_[ch];

				g = std::min(g, mg);
				g = std::max(g, -mg);

				drow[x] = ZP::to_unit(g);
			}
			srow += bpp;
			srow_ += bpp;
		}
	}
}


static void mask_to_byte(Image &maskbuf,PImage &mask,const PImage &_mask)
{
	if (n_channel(_mask) > 1 || depth(_mask) != 8)
	{
		maskbuf.construct(size(_mask), 1, 8);
		mask = maskbuf;

		switch (depth(_mask))
		{
		case 16:
			copy_channel<unsigned char, unsigned short>(
				mask, 0, _mask, 0, [](unsigned short v) {return v != 0 ? 255 : 0; });
			break;
		case 32:
			copy_channel<unsigned char, unsigned int>(
				mask, 0, _mask, 0, [](unsigned int v) {return v != 0 ? 255 : 0; });
			break;
		default:
			copy_channel<unsigned char, unsigned char>(
				mask, 0, _mask, 0, [](unsigned char v) {return v != 0 ? 255 : 0; });
		}
	}
	else
	{
		mask = _mask;
	}

}

void mask_merg(
	const bear::PImage &dst,
	const bear::PImage &src1,
	const bear::PImage &src2,
	const bear::PImage &mask)
{
	if (3 == src1.n_channel)
	{
		for (int y = 0; y < height(dst); ++y)
		{
			unsigned char * drow = (unsigned char *)scanline(dst, y);
			unsigned char * row1 = (unsigned char *)scanline(src1, y);
			unsigned char * row2 = (unsigned char *)scanline(src2, y);
			unsigned char * mrow = (unsigned char *)scanline(mask, y);

			for (int x = 0; x < width(dst); ++x)
			{
				bool f = 0 != mrow[x];

				if (f)
				{
					drow[0] = row2[0];
					drow[1] = row2[1];
					drow[2] = row2[2];
				}
				else
				{
					drow[0] = row1[0];
					drow[1] = row1[1];
					drow[2] = row1[2];
				}

				drow += 3;
				row1 += 3;
				row2 += 3;
			}
		}
	}
	else
	{
		for (int y = 0; y < height(dst); ++y)
		{
			unsigned int * drow = (unsigned int *)scanline(dst, y);
			unsigned int * row1 = (unsigned int *)scanline(src1, y);
			unsigned int * row2 = (unsigned int *)scanline(src2, y);
			unsigned char * mrow = (unsigned char *)scanline(mask, y);

			int bpp = src1.n_channel;

			for (int x = 0; x < width(dst); ++x)
			{
				bool f = 0 != mrow[x];

				if (f)
				{
					*drow = *row2;
				}
				else
				{
					*drow = *row1;
				}

				++drow;
				++row1;
				++row2;
			}
		}
	}
}

template<typename Unit>
static void stiching(const PImage &dst, const PImage ds, unsigned int ch)
{
	for (int y = 0; y < height(dst); ++y)
	{
		unsigned char * drow = (unsigned char *)scanline(dst, y);
		Unit * srow = (Unit *)scanline(ds, y);

		for (int x = 0; x < width(dst); ++x)
		{
			_ZERO_PS<Unit>::to_char(drow[ch], srow[x]);

			drow += dst.n_channel;
		}
	}
}

template<typename Unit>
static void _poisson_stiching_merged(
	const bear::PImage &dst,
	const bear::PImage &src,
	const bear::PImage &mask,
	unsigned int format,
	PStichingParam param)
{
	Image dx(size(dst), 1, sizeof(Unit) << 3);
	Image dy(size(dst), 1, sizeof(Unit) << 3);
	Image ds(size(dst), 1, sizeof(Unit) << 3);

	int chl[3];

	chl[0] = FI_RED(format);
	chl[1] = FI_GREEN(format);
	chl[2] = FI_BLUE(format);

	for (int i = 0; i < 3; ++i)
	{
		x_d<Unit>(dx, src, chl[i], mask, param.max_grandient);
		y_d<Unit>(dy, src, chl[i], mask, param.max_grandient);

		dxy_poisson_solver(ds, dx, dy, param.iteration_time, param.base_level);

		stiching<Unit>(dst, ds, chl[i]);
	}
}

void poisson_stiching_merged(
	const bear::PImage &dst,
	const bear::PImage &src,
	const bear::PImage &_mask,
	unsigned int format,
	PStichingParam param)
{
	Image maskbuf;
	PImage mask;
	mask_to_byte(maskbuf, mask, _mask);

	if (param.float_precision)
	{
		_poisson_stiching_merged<float>(dst, src, mask, format, param);
	}
	else
	{
		_poisson_stiching_merged<unsigned short>(dst, src, mask, format, param);
	}
}


void poisson_stiching(
	const bear::PImage &dst,
	const bear::PImage &src1,
	const bear::PImage &src2,
	const bear::PImage &_mask,
	unsigned int format,
	PStichingParam param)
{
	Image maskbuf;
	PImage mask;
	mask_to_byte(maskbuf, mask, _mask);

	mask_merg(src1, src1, src2, mask);

	poisson_stiching_merged(dst, src1, mask, format, param);
}