#include "possion_stiching.h"
#include "possion_solver.h"

#include <assert.h>
#include <cmath>

#include <algorithm>
#include <numeric>
#include <vector>

#include "utility.hpp"

using namespace bear;
using namespace std;


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
	int bpp = src.n_channel;

	fill(dx, (Unit)ZP::run());

	for (int y = 0; y < height(src); ++y)
	{
		Unit * drow = (Unit *)scanline(dx, y);
		unsigned char * srow = (unsigned char *)scanline(src, y);
		unsigned char * mrow = (unsigned char *)scanline(mask, y);

		srow += bpp;

		for (int x = 1; x < width(src); ++x)
		{
			bool f1 = 0 != mrow[x - 1];
			bool f2 = 0 != mrow[x];

			if (f1 != f2)
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
	int bpp = src.n_channel;

	fill(dy, (Unit)ZP::run());

	for (int y = 1; y < height(src); ++y)
	{
		Unit * drow = (Unit *)scanline(dy, y);
		unsigned char * srow = (unsigned char *)scanline(src, y);
		unsigned char * mrow = (unsigned char *)scanline(mask, y);

		unsigned char * srow_ = (unsigned char *)scanline(src, y - 1);
		unsigned char * mrow_ = (unsigned char *)scanline(mask, y - 1);

		for (int x = 0; x < width(src); ++x)
		{
			bool f1 = 0 != mrow_[x];
			bool f2 = 0 != mrow[x];

			if (f1 != f2)
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

static void conv_hst(vector<unsigned int> &hst)
{
	unsigned int buf[5] = {
		0,0,0,hst[0],hst[1],
	};

	for (unsigned int i = 0; i < hst.size(); ++i)
	{
		move(buf + 1, buf + 5, buf);
		buf[4] = (i < hst.size() - 2) ? hst[i + 2] : 0;

		hst[i] = accumulate(buf, buf + 5, 0);
	}
}


struct LimiteBound
{
	unsigned int tt;
	int maxg, ming;
	vector<unsigned int> hst,convhst;

	void init()
	{
		hst.resize(511);
	}

	void clear()
	{
		tt = 0;
		memset(&hst[0], 0, hst.size() * sizeof(unsigned int));
	}

	void stat(int ss1, int ss0)
	{
		int lm = maxg >> 1;

		if (ss1 > lm || ss0 > lm)
		{
			unsigned int g = ss1 - ss0 + 255;
			++hst[g];
			++tt;
		}

	}

	void get_limite(
		float ignore,
		float relax)
	{
		unsigned int num = (unsigned int)(tt * (1.0f - ignore));
		if (num < 5)return;

		convhst = hst;

		for (int i = 0; i < 3; ++i)conv_hst(convhst);

		int me = (int)(max_element(convhst.begin(), convhst.end()) - convhst.begin());

		unsigned int acc = hst[me];

		int min_l = 0;
		int max_l = 0;

		for (int i = 1; i < 511; ++i)
		{
			if (acc >= num)break;

			if (me - i >= 0 && hst[me - i])
			{
				acc += hst[me - i];
				min_l = -i;
			}
			if (me + i < (int)hst.size() && hst[me + i])
			{
				acc += hst[me + i];
				max_l = i;
			}
		}

		min_l = (int)(min_l*relax);
		max_l = (int)(max_l*relax);

		ming = max(ming, me + min_l - 255);
		maxg = min(maxg, me + max_l - 255);

		//ming = min(ming, 0);
		//maxg = max(maxg, 0);
	}
};

template<typename Unit>
static void x_d(
	const PImage &dx,
	const PImage &src,
	unsigned int ch,
	MakeMask &mm,
	int max_grandient)
{
	typedef _ZERO_PS<Unit> ZP;

	int bpp = src.n_channel;

	fill(dx, (Unit)ZP::run());

	LimiteBound lb;

	if (mm.fac.enable)lb.init();

	for (unsigned int y_b = 0; y_b < mm.h_border.size(); ++y_b)
	{
		int start_y = y_b ? mm.h_border[y_b - 1] : 0;
		int end_y = mm.h_border[y_b];

		int ph = end_y - start_y;
		if (!ph)continue;

		for (unsigned int x_b = 0; x_b < mm.w_border.size() - 1; ++x_b)
		{
			int x = mm.w_border[x_b];

			lb.maxg = max_grandient & 255;
			lb.ming = -lb.maxg;

			if (mm.fac.enable)
			{
				lb.clear();

				unsigned char * start_src = pick_pixel(src, x - 1, start_y);

				for (int i = 0; i < ph; ++i)
				{
					lb.stat(start_src[bpp + ch], start_src[ch]);
					start_src += src.width_step;
				}

				lb.get_limite(mm.fac.ignore, mm.fac.relax);
			}

			unsigned char * start_src = pick_pixel(src, x - 1, start_y);
			unsigned char * start_dst = pick_pixel(dx, x, start_y);

			for (int i = 0; i < ph; ++i)
			{
				int g = start_src[bpp + ch] - start_src[ch];

				g = std::min(g, lb.maxg);
				g = std::max(g, lb.ming);

				*((Unit *)start_dst) = ZP::to_unit(g);

				start_src += src.width_step;
				start_dst += dx.width_step;
			}
		}
	}
}

template<typename Unit>
static void y_d(
	const PImage &dy,
	const PImage &src,
	unsigned int ch,
	MakeMask &mm,
	int max_grandient)
{
	typedef _ZERO_PS<Unit> ZP;

	int bpp = src.n_channel;

	fill(dy, (Unit)ZP::run());

	LimiteBound lb;

	if (mm.fac.enable)lb.init();

	for (unsigned int x_b = 0; x_b < mm.w_border.size(); ++x_b)
	{
		int start_x = x_b ? mm.w_border[x_b - 1] : 0;
		int end_x = mm.w_border[x_b];

		int ph = end_x - start_x;
		if (!ph)continue;

		for (unsigned int y_b = 0; y_b < mm.h_border.size() - 1; ++y_b)
		{
			int y = mm.h_border[y_b];

			lb.maxg = max_grandient & 255;
			lb.ming = -lb.maxg;

			if (mm.fac.enable)
			{
				lb.clear();

				unsigned char * _start_src = pick_pixel(src, start_x, y - 1);
				unsigned char * start_src = pick_pixel(src, start_x, y);

				for (int i = 0; i < ph; ++i)
				{
					lb.stat(start_src[ch], _start_src[ch]);
					start_src += bpp;
					_start_src += bpp;
				}

				lb.get_limite(mm.fac.ignore, mm.fac.relax);
			}

			unsigned char * _start_src = pick_pixel(src, start_x, y - 1);
			unsigned char * start_src = pick_pixel(src, start_x, y);
			Unit * start_dst = (Unit *)pick_pixel(dy, start_x, y);

			for (int i = 0; i < ph; ++i)
			{
				int g = start_src[ch] - _start_src[ch];

				g = std::min(g, lb.maxg);
				g = std::max(g, lb.ming);

				*start_dst = ZP::to_unit(g);

				start_src += bpp;
				_start_src += bpp;
				++start_dst;
			}
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
		case 8:
			transform_channel<unsigned char, unsigned char>(
				mask, 0, _mask, 0, [](unsigned char v) {return v != 0 ? 255 : 0; });
			break;
		case 16:
			transform_channel<unsigned char, unsigned short>(
				mask, 0, _mask, 0, [](unsigned short v) {return v != 0 ? 255 : 0; });
			break;
		case 32:
			transform_channel<unsigned char, unsigned int>(
				mask, 0, _mask, 0, [](unsigned int v) {return v != 0 ? 255 : 0; });
			break;
		default:
			assert(false);
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
			//unsigned char u = 128;
			//_ZERO_PS<Unit>::to_char(u, srow[x]);
			//drow[ch] = u;

			_ZERO_PS<Unit>::to_char(drow[ch], srow[x]);

			drow += dst.n_channel;
		}
	}
}


inline void validate_border(vector<unsigned int> &border, unsigned int bound)
{
	border.erase(remove_if(border.begin(), border.end(), [bound](unsigned int b)
	{
		return b <= 0 || b >= bound;
	}), border.end());
}

template<typename Unit>
static void _poisson_stiching_merged(
	const bear::PImage &dst,
	const bear::PImage &src,
	MakeMask &mm,
	unsigned int format,
	PStichingParam param)
{
	Image maskbuf;
	PImage mask;

	if (mm.mask)
	{
		assert(
			width(src) == width(mm.mask) &&
			height(src) == height(mm.mask));
		mask_to_byte(maskbuf, mask, mm.mask);
	}
	else
	{
		validate_border(mm.w_border, width(src));
		validate_border(mm.h_border, height(src));

		sort(mm.w_border.begin(), mm.w_border.end());
		sort(mm.h_border.begin(), mm.h_border.end());

		mm.w_border.push_back(width(src));
		mm.h_border.push_back(height(src));
	}


	Image dx(size(dst), 1, sizeof(Unit) << 3);
	Image dy(size(dst), 1, sizeof(Unit) << 3);
	Image ds(size(dst), 1, sizeof(Unit) << 3);

	int chl[3];

	chl[0] = FI_RED(format);
	chl[1] = FI_GREEN(format);
	chl[2] = FI_BLUE(format);

	for (int i = 0; i < 3; ++i)
	{
		if (mm.mask)
		{
			x_d<Unit>(dx, src, chl[i], mask, param.max_grandient);
			y_d<Unit>(dy, src, chl[i], mask, param.max_grandient);
		}
		else
		{
			x_d<Unit>(dx, src, chl[i], mm, param.max_grandient);
			y_d<Unit>(dy, src, chl[i], mm, param.max_grandient);
		}


		dxy_poisson_solver(ds, dx, dy, param.iteration_time, param.base_level);

		stiching<Unit>(dst, ds, chl[i]);
	}
}

void poisson_stiching_merged(
	const bear::PImage &dst,
	const bear::PImage &src,
	MakeMask &&mask,
	unsigned int format,
	PStichingParam param)
{
	assert(
		n_channel(dst) == FI_BPP(format) &&
		n_channel(src) == FI_BPP(format) &&
		8 == depth(dst) &&
		8 == depth(src) &&
		width(src) == width(dst) &&
		height(src) == height(dst));

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