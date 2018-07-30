#include "possion_solver.h"
#include "utility.hpp"
#include "filter_convolution.hpp"

#include <assert.h>
#include <cmath>
#include <algorithm>

using namespace bear;


struct DownKernel
{
	enum { Size = 5 };
	enum { wi1 = 2048, wi2 = 8192, wi3 = 8192 + 4096, ofs = 16384, sft = 15 };
	constexpr static float w1()
	{
		return 0.0625f;
	}
	constexpr static float w2()
	{
		return 0.25f;
	}
	constexpr static float w3()
	{
		return 0.375f;
	}
};

struct UpKernel
{
	enum { Size = 5 };
	enum { wi1 = 4096, wi2 = 16384, wi3 = 16384 + 8192, ofs = 16384, sft = 15 };
	constexpr static float w1()
	{
		return 0.125f;
	}
	constexpr static float w2()
	{
		return 0.5f;
	}
	constexpr static float w3()
	{
		return 0.75f;
	}
};


static void set_zero_sign_u16(const PImage &img)
{
	for (int y = 0; y<img.height; ++y)
	{
		unsigned short * row = (unsigned short *)scanline(img, y);
		for (int x = 0; x<img.width; ++x)
		{
			row[x] = 32768;
		}
	}
}

static void iteration_u16(
	const PImage &dst,
	const PImage &dx,
	const PImage &dy,
	int iteration_time)
{
	for (int k = 0; k<iteration_time; ++k)
	{
		int y_flag = k & 1;

		for (int y = 0; y<dst.height; ++y)
		{
			int x_flag = (y + y_flag) & 1;

			unsigned short * xrow1 = (unsigned short *)scanline(dx, y - 1);
			unsigned short * xrow2 = (unsigned short *)scanline(dx, y);
			unsigned short * xrow3 = (unsigned short *)scanline(dx, y + 1);

			unsigned short * yrow1 = (unsigned short *)scanline(dy, y - 1);
			unsigned short * yrow2 = (unsigned short *)scanline(dy, y);
			unsigned short * yrow3 = (unsigned short *)scanline(dy, y + 1);


			unsigned short * drow1 = (unsigned short *)scanline(dst, y - 1);
			unsigned short * drow2 = (unsigned short *)scanline(dst, y);
			unsigned short * drow3 = (unsigned short *)scanline(dst, y + 1);

			if (0 == y)
			{
				int x = x_flag;

				if (!x)
				{
					drow2[x] = limiteU16((int)(
						32768 - xrow2[x + 1] + 32768 - yrow3[x]
						+ drow2[x + 1] + drow3[x] + 1
						) >> 1);

					x += 2;
				}
				int mx = dst.width - 1;
				for (; x<mx; x += 2)
				{
					drow2[x] = limiteU16((int)(
						xrow2[x] - xrow2[x + 1] + 32768 - yrow3[x]
						+ drow2[x - 1] + drow2[x + 1] + drow3[x] + 2
						) / 3);
				}
				if (x == mx)
				{
					drow2[x] = limiteU16((int)(
						xrow2[x] - 32768 + 32768 - yrow3[x]
						+ drow2[x - 1] + drow3[x] + 1
						) >> 1);
				}
			}
			else if(dst.height - 1 == y)
			{
				int x = x_flag;

				if (!x)
				{
					drow2[x] = limiteU16((int)(
						32768 - xrow2[x + 1] + yrow2[x] - 32768
						+ drow2[x + 1] + drow1[x] + 1
						) >> 1);

					x += 2;
				}
				int mx = dst.width - 1;
				for (; x<mx; x += 2)
				{
					drow2[x] = limiteU16((int)(
						xrow2[x] - xrow2[x + 1] + yrow2[x] - 32768
						+ drow2[x - 1] + drow2[x + 1] + drow1[x] + 2
						) / 3);
				}
				if (x == mx)
				{
					drow2[x] = limiteU16((int)(
						xrow2[x] - 32768 + yrow2[x] - 32768
						+ drow2[x - 1] + drow1[x] + 1
						) >> 1);
				}
			}
			else
			{
				int x = x_flag;

				if (!x)
				{
					drow2[x] = limiteU16((int)(
						32768 - xrow2[x + 1] + yrow2[x] - yrow3[x]
						+ drow2[x + 1] + drow1[x] + drow3[x] + 2
						) / 3);

					x += 2;
				}
				int mx = dst.width - 1;
				for (; x<mx; x += 2)
				{
					drow2[x] = limiteU16((int)(
						xrow2[x] - xrow2[x + 1] + yrow2[x] - yrow3[x]
						+ drow2[x - 1] + drow2[x + 1] + drow1[x] + drow3[x] + 2
						) >> 2);
				}
				if (x == mx)
				{
					drow2[x] = limiteU16((int)(
						xrow2[x] - 32768 + yrow2[x] - yrow3[x]
						+ drow2[x - 1] + drow1[x] + drow3[x] + 2
						) / 3);
				}
			}
		}
	}
}

static void scale_recursion_u16(
	const PImage & dst,
	const PImage &dx,
	const PImage &dy,
	int iteration_time, unsigned int n_layer)
{
	if (!n_layer)
	{
		set_zero_sign_u16(dst);
		return;
	}

	PSize subsz = size_down<5>(dst.size());

	Image sub_dx(subsz, 1, 16);
	Image sub_dy(subsz, 1, 16);
	Image sub_dst(subsz, 1, 16);

	DownFilter<DownKernel>::template run_dx<unsigned short>(sub_dx, dx);
	DownFilter<DownKernel>::template run_dy<unsigned short>(sub_dy, dy);


	scale_recursion_u16(sub_dst, sub_dx,sub_dy, iteration_time, n_layer - 1);

	UpFilter<UpKernel>::template run<unsigned short>(dst,sub_dst);

	iteration_u16(dst, dx, dy, iteration_time);
}

void dxy_poisson_solver_u16(
	const bear::PImage &_dst,
	const bear::PImage &dx,
	const bear::PImage &dy,
	unsigned int iteration_time,
	int base_level)
{
	Image dbuf;
	PImage dst;
	if (scanline(dst, 0) == scanline(dx, 0) || scanline(dst, 0) == scanline(dy, 0))
	{
		dbuf.construct(_dst.size(), 1, 16);
		dst = dbuf;
	}
	else
	{
		dst = _dst;
	}

	float os = sqrt((float)dst.width * dst.width + dst.height * dst.height);


	int level = 0;

	while (os > 1.0f)
	{
		++level;
		os *= 0.5f;
	}

	level -= base_level;

	level = std::max(1, level);

	scale_recursion_u16(dst, dx, dy, iteration_time, level);
}

void dxy_poisson_solver(
	const bear::PImage &dst,
	const bear::PImage &dx,
	const bear::PImage &dy,
	unsigned int iteration_time,
	int base_level)
{
	if (dst.depth == 16 && dx.depth == 16 && dy.depth == 16)
	{
		dxy_poisson_solver_u16(dst, dx, dy, iteration_time, base_level);

		return;
	}

	assert(false);
}