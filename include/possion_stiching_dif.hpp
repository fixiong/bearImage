#ifndef _POSSION_STICHING_DIF_HPP
#define _POSSION_STICHING_DIF_HPP

#include "possion_stiching.h"
#include <memory.h>


using namespace bear;
using namespace std;


static void x_d(const image_ptr<unsigned short,1> &img)
{
	for (size_t y = 0; y < img.height(); ++y)
	{
		auto row = img[y];
		for (size_t x = img.width() - 1; x > 0; --x)
		{
			row[x] = row[x] - row[x - 1] + 32768;
		}

		row[0] = 32768;
	}
}


static void y_d(const image_ptr<unsigned short, 1> &img)
{
	for (size_t y = img.height() - 1; y >= 0; --y)
	{
		auto row_1 = img[y - 1];
		auto row = img[y];
		for (size_t x = 0; x < img.width(); ++x)
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

	static float from_unit(unsigned char v1,unsigned char v2)
	{
		return (float)((int)v2 - (int)v1);
	}

	static void to_unit(unsigned char &dst, float v)
	{
		dst = limiteU8((int)(dst - v + 0.5f));
	}

	static float from_unit(unsigned short v1, unsigned short v2)
	{
		return (float)((int)v2 - (int)v1);
	}

	static void to_unit(unsigned short &dst, float v)
	{
		dst = limiteU16((int)(dst - v + 0.5f));
	}
};

template<>
struct _ZERO_PS<unsigned short>
{
	constexpr static unsigned short run()
	{
		return 32768;
	}

	static unsigned short from_unit(unsigned char v1, unsigned char v2)
	{
		return (256 + v2 - v1) << 7;
	}

	static void to_unit(unsigned char &dst, unsigned short v)
	{
		dst = limiteU8(dst - round_shift<7>((int)v - 32768));
	}

	static unsigned short from_unit(unsigned short v1, unsigned short v2)
	{
		return (65536 + v2 - v1) >> 1;
	}

	static void to_unit(unsigned short &dst, unsigned short v)
	{
		dst = limiteU8(dst - ((int)(v << 1) - 65535));
	}
};


template<typename DstUnit, typename SrcUnit, typename ZP>
static void x_d(
	image_ptr<DstUnit,1> dx,
	tensor_ptr<SrcUnit,3> src,
	unsigned int ch,
	image_ptr<unsigned char,1> mask,
	ZP && zp)
{
	dx.fill(zp.run());

	for (size_t y = 0; y < size_at<0>(src); ++y)
	{
		auto drow = dx[y];
		auto srow = src[y];
		auto mrow = mask[y];

		for (size_t x = 1; x < size_at<1>(src); ++x)
		{
			bool f1 = 0 != mrow[x - 1];
			bool f2 = 0 != mrow[x];

			if (f1 != f2)
			{
				drow[x] = zp.from_unit(srow[x - 1][ch], srow[x][ch]);
			}
		}
	}
}

template<typename DstUnit, typename SrcUnit, typename ZP>
static void y_d(
	image_ptr<DstUnit, 1> &dy,
	tensor_ptr<SrcUnit, 3> &src,
	unsigned int ch,
	image_ptr<unsigned char, 1> mask,
	ZP && zp)
{
	dy.fill(zp.run());

	for (int y = 1; y < size_at<0>(src); ++y)
	{
		auto drow = dy[y];
		auto srow = src[y];
		auto mrow = mask[y];

		auto srow_ = src[y - 1];
		auto mrow_ = mask[y - 1];

		for (int x = 0; x < size_at<1>(src); ++x)
		{
			bool f1 = 0 != mrow_[x];
			bool f2 = 0 != mrow[x];

			if (f1 != f2)
			{
				drow[x] = zp.from_unit(srow_[x],srow[x]);
			}
		}
	}
}

template<typename DstUnit, typename SrcUnit, typename ZP>
static void x_d_p(
	image_ptr<DstUnit,1> dx,
	tensor_ptr<SrcUnit, 3> src1,
	tensor_ptr<SrcUnit, 3> src2,
	unsigned int ch,
	ZP && zp)
{
	unsigned int rd = dx.width() / 2;
	double bs = 1.0 / dx.width();

	for (size_t y = 0; y < dx.height(); ++y)
	{
		double t = 0.0;
		for (size_t x = 0; x < dx.width(); ++x)
		{
			t += zp.from_unit(src1[y][x][ch], src2[y][x][ch]);
		}

		dx[y][rd] = (DstUnit)floor(t * bs + 0.5);
	}
}

template<typename DstUnit, typename SrcUnit, typename ZP>
static void y_d_p(
	image_ptr<DstUnit,1> dy,
	tensor_ptr<SrcUnit, 3> src1,
	tensor_ptr<SrcUnit, 3> src2,
	unsigned int ch,
	ZP && zp)
{
	unsigned int rd = dy.height() / 2;
	double bs = 1.0 / dx.height();

	for (size_t x = 0; x < dy.width(); ++x)
	{
		double t = 0.0;
		for (size_t y = 0; y < dy.height(); ++y)
		{
			t += zp.from_unit(src1[y][x][ch], src2[y][x][ch]);
		}

		dx[rd][x] = (DstUnit)floor(t * bs + 0.5);
	}
}


#endif