#include "../include/possion_solver.h"
#include "../include/utility.h"
#include "filter_convolution.hpp"

#include "../../bear/include/ptr_algorism.h"

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

template<typename Unit>
struct Iteration {};


template<>
struct Iteration<unsigned short>
{
	static void run(
		image_ptr<unsigned short, 1> dst,
		image_ptr<unsigned short, 1> dx,
		image_ptr<unsigned short, 1> dy,
		int iteration_time);
};

void Iteration<unsigned short>::run(
	image_ptr<unsigned short, 1> dst,
	image_ptr<unsigned short, 1> dx,
	image_ptr<unsigned short, 1> dy,
	int iteration_time)
{
	for (int k = 0; k<iteration_time; ++k)
	{
		int y_flag = k & 1;

		for (size_t y = 0; y<dst.height(); ++y)
		{
			int x_flag = (y + y_flag) & 1;

			auto ym = y - 1 < dst.height() ? y - 1 : 0;
			auto yp = y + 1 < dst.height() ? y + 1 : dst.height() - 1;

			auto xrow1 = dx[ym];
			auto xrow2 = dx[y];
			auto xrow3 = dx[yp];

			auto yrow1 = dy[ym];
			auto yrow2 = dy[y];
			auto yrow3 = dy[yp];


			auto drow1 = dst[ym];
			auto drow2 = dst[y];
			auto drow3 = dst[yp];

			if (0 == y)
			{
				size_t x = x_flag;

				if (!x)
				{
					drow2[x] = limiteU16((int)(
						32768 - xrow2[x + 1] + 32768 - yrow3[x]
						+ drow2[x + 1] + drow3[x] + 1
						) >> 1);

					x += 2;
				}
				size_t mx = dst.width() - 1;
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
			else if(dst.height() - 1 == y)
			{
				size_t x = x_flag;

				if (!x)
				{
					drow2[x] = limiteU16((int)(
						32768 - xrow2[x + 1] + yrow2[x] - 32768
						+ drow2[x + 1] + drow1[x] + 1
						) >> 1);

					x += 2;
				}
				size_t mx = dst.width() - 1;
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
				size_t x = x_flag;

				if (!x)
				{
					drow2[x] = limiteU16((int)(
						32768 - xrow2[x + 1] + yrow2[x] - yrow3[x]
						+ drow2[x + 1] + drow1[x] + drow3[x] + 2
						) / 3);

					x += 2;
				}
				size_t mx = dst.width() - 1;
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


template<>
struct Iteration<float>
{
	static void run(
		image_ptr<float, 1> dst,
		image_ptr<float, 1> dx,
		image_ptr<float, 1> dy,
		int iteration_time);
};

void Iteration<float>::run(
	image_ptr<float, 1> dst,
	image_ptr<float, 1> dx,
	image_ptr<float, 1> dy,
	int iteration_time)
{
	image_ptr<float, 1> lp = dx;

	for (size_t y = 0; y<dst.height(); ++y)
	{
		auto yp = y + 1 < dst.height() ? y + 1 : dst.height() - 1;

		auto xrow = dx[y];

		auto yrow1 = dy[y];
		auto yrow2 = dy[yp];

		auto drow = lp[y];

		if (y == dst.height() - 1)
		{
			auto mx = dst.width() - 1;
			for (size_t x = 0; x<mx; ++x)
			{
				drow[x] = xrow[x] - xrow[x + 1] + yrow1[x];
			}
			drow[mx] = xrow[mx] + yrow1[mx];
		}
		else
		{
			auto mx = dst.width() - 1;
			for (size_t x = 0; x<mx; ++x)
			{
				drow[x] = xrow[x] - xrow[x + 1] + yrow1[x] - yrow2[x];
			}
			drow[mx] = xrow[mx] + yrow1[mx] - yrow2[mx];

		}
	}

	for (int k = 0; k<iteration_time; ++k)
	{
		int y_flag = k & 1;

		for (size_t y = 0; y<dst.height(); ++y)
		{
			int x_flag = (y + y_flag) & 1;

			auto ym = y - 1 < dst.height() ? y - 1 : 0;
			auto yp = y + 1 < dst.height() ? y + 1 : dst.height() - 1;

			auto lprow = lp[y];

			auto drow1 = dst[ym];
			auto drow2 = dst[y];
			auto drow3 = dst[yp];

			if (0 == y)
			{
				size_t x = x_flag;

				if (!x)
				{
					drow2[x] = (lprow[x] + drow2[x + 1] + drow3[x]) * 0.5f;

					x += 2;
				}
				auto mx = dst.width() - 1;
				for (; x<mx; x += 2)
				{
					drow2[x] = (lprow[x] +
						drow2[x - 1] + drow2[x + 1] + drow3[x]) * (1.0f / 3.0f);
				}
				if (x == mx)
				{
					drow2[x] = (lprow[x] + drow2[x - 1] + drow3[x]) * 0.5f;
				}
			}
			else if (dst.height() - 1 == y)
			{
				size_t x = x_flag;

				if (!x)
				{
					drow2[x] = (lprow[x] + drow2[x + 1] + drow1[x]) * 0.5f;

					x += 2;
				}
				auto mx = dst.width() - 1;
				for (; x<mx; x += 2)
				{
					drow2[x] = (lprow[x] +
						drow2[x - 1] + drow2[x + 1] + drow1[x]) * (1.0f / 3.0f);
				}
				if (x == mx)
				{
					drow2[x] = (lprow[x] + drow2[x - 1] + drow1[x]) * 0.5f;
				}
			}
			else
			{
				size_t x = x_flag;

				if (!x)
				{
					drow2[x] = (lprow[x] +
						drow2[x + 1] + drow1[x] + drow3[x]) * (1.0f / 3.0f);

					x += 2;
				}
				auto mx = dst.width() - 1;
				for (; x<mx; x += 2)
				{
					drow2[x] = (lprow[x] +
						drow2[x - 1] + drow2[x + 1] + drow1[x] + drow3[x]) * 0.25f;
				}
				if (x == mx)
				{
					drow2[x] = (lprow[x] +
						drow2[x - 1] + drow1[x] + drow3[x]) * (1.0f / 3.0f);
				}
			}
		}
	}
}

template<typename Unit>
static void scale_recursion(
	image_ptr<Unit, 1> dst,
	image_ptr<Unit, 1> dx,
	image_ptr<Unit, 1> dy,
	int iteration_time, unsigned int n_layer)
{
	if (!n_layer)
	{
		dst.fill(_ZERO_UNIT<Unit>::run());
		return;
	}

	auto subsz = size_down<5>(size(dst));

	image<Unit, 1> sub_dx(subsz);
	image<Unit, 1> sub_dy(subsz);
	image<Unit, 1> sub_dst(subsz);

	DownFilter<DownKernel>::run_dx(to_ptr(sub_dx), dx);
	DownFilter<DownKernel>::run_dy(to_ptr(sub_dy), dy);


	scale_recursion<Unit>(sub_dst, sub_dx,sub_dy, iteration_time, n_layer - 1);

	UpFilter<UpKernel>::run(dst,to_ptr(sub_dst));

	Iteration<Unit>::run(dst, dx, dy, iteration_time);
}


template<typename Unit>
static void scale_recursion(
	image_ptr<Unit, 1> dst,
	image_ptr<Unit, 1> dx,
	image_ptr<Unit, 1> dy,
	image_ptr<Unit, 1> xborder,
	image_ptr<Unit, 1> yborder,
	int iteration_time, unsigned int n_layer)
{
	if (!n_layer)
	{
		dst.fill(_ZERO_UNIT<Unit>::run());
		return;
	}

	auto subsz = size_down<5>(size(dst));

	image<Unit, 1> sub_dx(subsz);
	image<Unit, 1> sub_dy(subsz);
	image<Unit, 1> sub_dst(subsz);
	image<Unit, 1> sub_xborder(2, subsz.height);
	image<Unit, 1> sub_yborder(subsz.width, 2);

	DownFilter<DownKernel>::run_dx(to_ptr(sub_dx), dx);
	DownFilter<DownKernel>::run_dy(to_ptr(sub_dy), dy);


	DownFilter<DownKernel>::run_bx(to_ptr(sub_xborder), xborder);
	DownFilter<DownKernel>::run_by(to_ptr(sub_yborder), yborder);


	scale_recursion<Unit>(sub_dst, sub_dx, sub_dy, sub_xborder, sub_yborder, iteration_time, n_layer - 1);

	UpFilter<UpKernel>::run(dst, to_ptr(sub_dst));

	Iteration<Unit>::run(dst, dx, dy, xborder, yborder, iteration_time);
}

template<typename Unit>
void dxy_poisson_solver_inner(
	image_ptr<Unit, 1> _dst,
	image_ptr<Unit, 1> dx,
	image_ptr<Unit, 1> dy,
	image_ptr<Unit, 1> xborder,
	image_ptr<Unit, 1> yborder,
	unsigned int iteration_time,
	int base_level)
{
	//copy(_dst, dy);
	//return;

	image<Unit,1> dbuf;
	image_ptr<Unit,1> dst;
	if (_dst.data() == dx.data() || _dst.data() == dy.data())
	{
		dbuf = image<Unit, 1>(size(_dst));
		dst = dbuf;
	}
	else
	{
		dst = _dst;
	}

	float os = sqrt((float)dst.width() * dst.width() + dst.height() * dst.height());


	int level = 0;

	while (os > 1.0f)
	{
		++level;
		os *= 0.5f;
	}

	level -= base_level;

	level = std::max(1, level);

	scale_recursion(dst, dx, dy, iteration_time, level);

	if (dbuf.size())copy(_dst, dst);
}

void dxy_poisson_solver(
	dynamic_image_ptr dst,
	dynamic_image_ptr dx,
	dynamic_image_ptr dy,
	unsigned int iteration_time,
	int base_level)
{

	if (dst.elm_size() == 2)
	{
		dxy_poisson_solver_inner(
			image_ptr<unsigned short, 1>(dst),
			image_ptr<unsigned short, 1>(dx),
			image_ptr<unsigned short, 1>(dy),
			image_ptr<unsigned short, 1>(),
			image_ptr<unsigned short, 1>(),
			iteration_time, base_level);
	}
	else if (dst.elm_size() == 4)
	{
		dxy_poisson_solver_inner(
			image_ptr<float, 1>(dst),
			image_ptr<float, 1>(dx),
			image_ptr<float, 1>(dy),
			image_ptr<float, 1>(),
			image_ptr<float, 1>(),
			iteration_time, base_level);
	}
	else
	{
		throw bear_exception(exception_type::other_error, "unsupported image type!");
	}
}


void dxy_poisson_solver(
	bear::dynamic_image_ptr dst,
	bear::dynamic_image_ptr dx,
	bear::dynamic_image_ptr dy,
	bear::dynamic_image_ptr xborder,
	bear::dynamic_image_ptr yborder,
	unsigned int iteration_time,
	int base_level)
{
	if (dst.elm_size() == 2)
	{
		dxy_poisson_solver_inner(
			image_ptr<unsigned short, 1>(dst),
			image_ptr<unsigned short, 1>(dx),
			image_ptr<unsigned short, 1>(dy),
			image_ptr<unsigned short, 1>(xborder),
			image_ptr<unsigned short, 1>(yborder),
			iteration_time, base_level);
	}
	else if (dst.elm_size() == 4)
	{
		dxy_poisson_solver_inner(
			image_ptr<float, 1>(dst),
			image_ptr<float, 1>(dx),
			image_ptr<float, 1>(dy),
			image_ptr<float, 1>(xborder),
			image_ptr<float, 1>(yborder),
			iteration_time, base_level);
	}
	else
	{
		throw bear_exception(exception_type::other_error, "unsupported image type!");
	}
}