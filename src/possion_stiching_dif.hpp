#ifndef _POSSION_STICHING_DIF_HPP
#define _POSSION_STICHING_DIF_HPP

#include <bear/image.h>


using namespace bear;
using namespace std;



template<typename T>
struct _ZERO_PS
{
};

template<>
struct _ZERO_PS<float>
{
	constexpr static float run()
	{
		return 0.0f;
	}

	constexpr static float one()
	{
		return 1.0f;
	}

	static float limite_zero(float v)
	{
		return v;
	}

	static float limite(float v)
	{
		return v;
	}

	static float minus(float v)
	{
		return -v;
	}

	static float acc(float v1, float v2)
	{
		return v1 + v2;
	}

	static float from_unit(unsigned char v1, unsigned char v2, unsigned char rv1, unsigned char rv2)
	{
		return (float)((int)v2 - (int)v1 - ((int)rv2 - (int)rv1));
	}

	static void to_unit(unsigned char& dst, float v)
	{
		dst = limiteU8((int)(dst - v + 0.5f));
	}

	static float from_unit(unsigned short v1, unsigned short v2, unsigned short rv1, unsigned short rv2)
	{
		return (float)((int)v2 - (int)v1 - ((int)rv2 - (int)rv1));
	}

	static void to_unit(unsigned short& dst, float v)
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

	constexpr static unsigned one()
	{
		return 65535;
	}

	static unsigned short limite_zero(int v)
	{
		return limiteU16(v + 32768);
	}

	static unsigned short limite(int v)
	{
		return limiteU16(v);
	}

	static unsigned short minus(unsigned short v)
	{
		return 65536 - v;
	}

	static int acc(int v1, int v2)
	{
		return v1 + v2 - 32768;
	}

	static unsigned short from_unit(unsigned char v1, unsigned char v2, unsigned char rv1, unsigned char rv2)
	{
		return (256 + v2 - v1 - rv2 + rv1) << 7;
	}

	static void to_unit(unsigned char& dst, unsigned short v)
	{
		dst = limiteU8(dst - round_shift<7>((int)v - 32768));
	}

	static unsigned short from_unit(unsigned short v1, unsigned short v2, unsigned short rv1, unsigned short rv2)
	{
		return (65536 + v2 - v1 - rv2 + rv1) >> 1;
	}

	static void to_unit(unsigned short& dst, unsigned short v)
	{
		dst = limiteU16(dst - ((int)(v << 1) - 65535));
	}
};


template<typename Dst, typename Src, typename Ref, typename ZP>
static void x_d(
	Dst dx,
	Src src,
	Ref ref,
	size_t ch,
	const_image_ptr<unsigned char, 1> mask,
	ZP&& zp)
{
	dx.fill(zp.run());

	for (size_t y = 0; y < height(src); ++y)
	{
		auto drow = dx[y];
		auto srow = src[y];
		auto rrow = ref[y];
		auto mrow = mask[y];

		for (size_t x = 1; x < width(src); ++x)
		{
			bool f1 = 0 != mrow[x - 1];
			bool f2 = 0 != mrow[x];

			if (f1 != f2)
			{
				drow[x] = std::forward<ZP>(zp).from_unit(srow[x - 1][ch], srow[x][ch], rrow[x - 1][ch], rrow[x][ch]);
			}
		}
	}
}

template<typename Dst, typename Src, typename Ref, typename ZP>
static void y_d(
	Dst dy,
	Src src,
	Ref ref,
	size_t ch,
	const_image_ptr<unsigned char, 1> mask,
	ZP&& zp)
{
	dy.fill(zp.run());

	for (size_t y = 1; y < height(src); ++y)
	{
		auto drow = dy[y];
		auto srow = src[y];
		auto rrow = ref[y];
		auto mrow = mask[y];

		auto srow_ = src[y - 1];
		auto rrow_ = ref[y - 1];
		auto mrow_ = mask[y - 1];

		for (size_t x = 0; x < width(src); ++x)
		{
			bool f1 = 0 != mrow_[x];
			bool f2 = 0 != mrow[x];

			if (f1 != f2)
			{
				drow[x] = std::forward<ZP>(zp).from_unit(srow_[x][ch], srow[x][ch], rrow_[x][ch], rrow[x][ch]);
			}
		}
	}
}

template<typename Dst, typename Src, typename ZP>
static void x_d_p(
	Dst dx,
	Src src1,
	Src src2,
	size_t ch,
	ZP&& zp)
{
	auto rd = dx.width() / 2;
	double bs = 1.0 / dx.width();

	for (size_t y = 0; y < dx.height(); ++y)
	{
		double t = 0.0;
		for (size_t x = 0; x < dx.width(); ++x)
		{
			t += std::forward<ZP>(zp).from_unit(src1[y][x][ch], src2[y][x][ch], 0, 0);
		}

		dx[y][rd] = (typename Dst::elm_type)floor(t * bs + 0.5);
	}
}

template<typename Dst, typename Src, typename ZP>
static void y_d_p(
	Dst dy,
	Src src1,
	Src src2,
	size_t ch,
	ZP&& zp)
{
	auto rd = dy.height() / 2;
	double bs = 1.0 / dy.height();

	for (size_t x = 0; x < dy.width(); ++x)
	{
		double t = 0.0;
		for (size_t y = 0; y < dy.height(); ++y)
		{
			t += std::forward<ZP>(zp).from_unit(src1[y][x][ch], src2[y][x][ch], 0, 0);
		}

		dy[rd][x] = (typename Dst::elm_type)floor(t * bs + 0.5);
	}
}


#endif