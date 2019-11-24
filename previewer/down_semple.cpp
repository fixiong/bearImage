#include "down_semple.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace bear;

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

template<typename D, typename ST, size_t SD, typename Fun>
inline void down_semple_line(D && dst, const array<ST, SD> &src, Fun && fun)
{
	for (auto i = 0; i < SD; ++i)
	{
		forward<Fun>(fun)(forward<D>(dst)[i], src[i]);
	}
}

inline array<size_t, 3> create_tmp(const array_ptr<array<size_t, 3>> &)
{
	return array<size_t, 3>();
}

inline vector<array<size_t, 3>> create_tmp(const image_ptr<unsigned char, 3> &dst)
{
	return vector<array<size_t, 3>>(width(dst));
}

template<typename D, typename S, typename Fun, typename ... F>
void down_semple_line(D && dst, S &&src, Fun && fun, float _left, float _fac , F ... other_fac)
{
	pos_t ileft = (pos_t)(_left * 256.0f + 0.5f);
	pos_t ifac = (pos_t)(65536.0f / _fac + 0.5f);
	pos_t istep = (ifac + 128) >> 8;

	auto dp = create_tmp(dst);

	for (size_t i = 0; i < dst.size(); ++i)
	{
		pos_t istart = ileft + ((i * ifac + 128) >> 8);
		pos_t index = istart >> 8;
		size_t start_weight = 256 - (istart & 255);

		pos_t iend = istart + istep;
		pos_t max_index = iend >> 8;
		size_t end_weight = 256 - (iend & 255);

		down_semple_line(dp, clip_at(src, index), [=](auto &d, auto s) {
			d = s * start_weight;
		}, other_fac ...);
		++index;

		while (index < max_index)
		{
			down_semple_line(dp, clip_at(src, index),[](auto &d, auto s) {
				d += s << 8;
			}, other_fac ...);
			++index;
		}

		down_semple_line(dp, clip_at(src, index),[=](auto &d, auto s) {
			d += (s << 8) - s * end_weight;
		}, other_fac ...);

		map_function(forward<Fun>(fun), forward<D>(dst)[i], to_ptr(dp));
	}
}

void down_semple(
	image_ptr<unsigned char, 3> dst,
	const_image_ptr<unsigned char, 3> src,
	float x_left,
	float x_fac,
	float y_left,
	float y_fac)
{
	float fac = (x_fac * y_fac * 65536);
	size_t sft = 0;

	while (fac < 32768.0f)
	{
		fac *= 2.0f;
		sft += 1;
	}
	unsigned long long ifac = (unsigned long long)(fac + 0.5f);

	down_semple_line(dst, src,[=](auto &dst, auto src) {
		dst = (unsigned char)((unsigned long long)(src * ifac) >> (32 + sft));
	}, y_left, y_fac, x_left, x_fac);
}