#include "down_semple.h"
#include <cmath>
#include <bear/ptr_numeric.h>
#include <iostream>

using namespace std;
using namespace bear;


using image_t = image<unsigned char, 3>;

template<typename T>
inline auto clip_at(T && p, pos_t index)
{
	if ((size_t)index >= p.size())
	{
		if (index < 0)
		{
			return forward<T>(p).at(0);
		}
		else
		{
			return forward<T>(p).at(p.size() - 1);
		}
	}
	return forward<T>(p).at((size_t)index);
}

struct dummy {
	dummy & operator[](size_t) {
		return *this;
	}
};

template<typename D1, typename D2, typename S, typename Fun>
inline void down_semple_pixel(D1 && dst1, D2 && dst2, S &&src, Fun && fun)
{
	forward<Fun>(fun)(forward<D1>(dst1)[0], forward<D2>(dst2)[0], forward<S>(src)[0]);
	forward<Fun>(fun)(forward<D1>(dst1)[1], forward<D2>(dst2)[1], forward<S>(src)[1]);
	forward<Fun>(fun)(forward<D1>(dst1)[2], forward<D2>(dst2)[2], forward<S>(src)[2]);

}

template<typename D1, typename D2, typename Fun>
void down_semple_line(D1 && dst1, D2 && dst2, const_array_ptr<array<unsigned char, 3>> src, float first_pos, float fac, Fun && fun)
{
	float step = 1.0f / fac;
	float start = first_pos - step / 2.0f;

	array<size_t, 3> dp;

	pos_t istart = pos_t(start * 256.0f);
	pos_t current = istart >> 8;
	size_t weight = size_t(((current + 1) << 8) - istart);
	down_semple_pixel(dp,dummy(), clip_at(src, current),[=](size_t &d, auto, unsigned char s) {
		d = s * weight;
	});
	++current;

	for (size_t x = 0; x < dst1.size(); ++x)
	{
		start += step;
		istart = pos_t(start * 256.0f);
		pos_t next = (istart >> 8) + 1;


		while (current < next - 1)
		{
			down_semple_pixel(dp, dummy(), clip_at(src, current),[](size_t &d, auto, unsigned char s) {
				d += s << 8;
			});
			++current;
		}

		array<size_t, 3> ndp;

		size_t weight = size_t(((current + 1) << 8) - istart);
		down_semple_pixel(dp, ndp, clip_at(src, current),[=](size_t &d, size_t &n, unsigned char s) {
			d += s << 8;
			n = s * weight;
			d -= n;
		});
		++current;

		down_semple_pixel(forward<D1>(dst1)[x], forward<D2>(dst2)[x], dp, forward<Fun>(fun));

		dp = ndp;
	}
}


void down_semple(image_ptr<unsigned char, 3> dst, const_image_ptr<unsigned char, 3> src, float x_first_Pos, float x_fac, float y_first_Pos, float y_fac)
{
	float fac = x_fac * y_fac / 256.0f / 256.0f;
	float step = 1.0f / y_fac;
	float start = y_first_Pos - step / 2.0f;

	tensor<array<size_t, 3>, 2> tmp(2, width(dst));

	auto dp = tmp[0];
	auto ndp = tmp[1];


	pos_t istart = pos_t(start * 256.0f);
	pos_t current = istart >> 8;
	size_t weight = size_t(((current + 1) << 8) - istart);
	down_semple_line(dp, dummy(), clip_at(src, current), x_first_Pos, x_fac, [=](size_t &d, auto, size_t s) {
		d = s * weight;
	});
	++current;

	for (size_t y = 0; y < dst.size(); ++y)
	{
		start += step;
		istart = pos_t(start * 256.0f);
		pos_t next = (istart >> 8) + 1;


		while (current < next - 1)
		{
			down_semple_line(dp, dummy(), clip_at(src, current), x_first_Pos, x_fac, [](size_t &d, auto, size_t s) {
				d += s << 8;
			});
			++current;
		}

		size_t weight = size_t(((current + 1) << 8) - istart);
		down_semple_line(dp, ndp, clip_at(src, current), x_first_Pos, x_fac,[=](size_t &d, size_t &n, size_t s) {
			d += s << 8;
			n = s * weight;
			d -= n;
		});
		++current;

		auto dy = dst[y];

		for (size_t x = 0; x < dp.size(); ++x)
		{
			dy[x][0] = (unsigned char)(dp[x][0] * fac + 0.5f);
			dy[x][1] = (unsigned char)(dp[x][1] * fac + 0.5f);
			dy[x][2] = (unsigned char)(dp[x][2] * fac + 0.5f);
		}

		swap(dp, ndp);
	}

}