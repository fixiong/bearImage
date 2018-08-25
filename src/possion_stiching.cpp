#include "possion_stiching.h"
#include "possion_solver.h"

#include <assert.h>
#include <cmath>

#include <algorithm>
#include <numeric>
#include <vector>

#include "utility.hpp"

#include "possion_stiching_dif.hpp"

using namespace bear;
using namespace std;


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
			//drow[ch] = 255 - u;

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


template<typename Unit,typename Dx, typename Dy>
static void _poisson_stiching_inner(
	const bear::PImage &dst,
	const bear::PImage &src,
	Dx && get_dx,
	Dy && get_dy,
	unsigned int format,
	PStichingParam param)
{
	copy(dst, src);

	Image dx(size(dst), 1, sizeof(Unit) << 3);
	Image dy(size(dst), 1, sizeof(Unit) << 3);
	Image ds(size(dst), 1, sizeof(Unit) << 3);

	int chl[3];

	chl[0] = FI_RED(format);
	chl[1] = FI_GREEN(format);
	chl[2] = FI_BLUE(format);

	for (int i = 0; i < 3; ++i)
	{
		std::forward<Dx>(get_dx)(dx, chl[i]);
		std::forward<Dy>(get_dy)(dy, chl[i]);

		dxy_poisson_solver(ds, dx, dy, param.iteration_time, param.base_level);

		stiching<Unit>(dst, ds, chl[i]);
	}
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

	_poisson_stiching_inner<Unit>(
		dst, src,
		[&](PImage dx, int ch)
	{
		if (mm.mask)
			x_d<Unit>(dx, src, ch, mask, param.max_grandient);
		else
			x_d<Unit>(dx, src, ch, mm, param.max_grandient);
	},
		[&](PImage dy, int ch)
	{
		if (mm.mask)
			y_d<Unit>(dy, src, ch, mask, param.max_grandient);
		else
			y_d<Unit>(dy, src, ch, mm, param.max_grandient);
	},
		format,
		param
	);
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


template<typename C>
static void for_each_img(const vector<vector<PImage>> &src,unsigned int rd, C && c)
{
	int nw = (int)src[0].size();
	int nh = (int)src.size();

	vector<int> bws;
	vector<int> bhs;

	int cpy = 0;
	for (int y = 0; y < nh; ++y)
	{
		assert(src[y].size() == nw);

		int rh = src[y][0].height;
		int bh = rh;
		int py = 0;

		if (0 != y)
		{
			bh -= rd;
			py += rd;
		}

		if (nh - 1 != y) bh -= rd;

		bhs.push_back(bh);

		int cpx = 0;
		for (int x = 0; x < nw; ++x)
		{
			assert(src[y][x].height == rh);

			int rw = src[y][x].width;
			int bw = rw;
			int px = 0;

			if (0 != x)
			{
				bw -= rd;
				px += rd;
			}

			if (nw - 1 != x) bw -= rd;

			if (y)assert(bw == bws[x]);
			else bws.push_back(bw);

			forward<C>(c)(
				src, x, y,
				PSize(bw, bh),
				PPoint(px, py),
				PPoint(cpx, cpy));

			cpx += bw;
		}

		cpy += bh;
	}

}

template<typename Unit>
static void _poisson_stiching_m(
	const bear::PImage &dst,
	const vector<vector<PImage>> &src,
	unsigned int rd,
	unsigned int format,
	PStichingParam param)
{

	_poisson_stiching_inner<Unit>(
		dst, dst,
		[&src, &param, rd](PImage dx, int ch)
	{
		clear_d<Unit>(dx);
		for_each_img(src, rd, [&dx, ch, &param, rd](
			const vector<vector<PImage>> &src,
			int x, int y,
			PSize bs,
			PPoint p,
			PPoint cp
			)
		{
			if (src[y].size() - 1 == x)return;

			auto img1 = src[y][x];
			auto img2 = src[y][x + 1];

			img1 = clip_image(img1, img1.width - rd * 2, p.y, rd * 2, bs.height);
			img2 = clip_image(img2, 0, p.y, rd * 2, bs.height);

			auto dst = clip_image(dx, cp.x + bs.width - rd, cp.y, rd * 2, bs.height);

			x_d_p<Unit>(dst, img1, img2, ch, param.max_grandient);

		});
	},
		[&src, &param, rd](PImage dy, int ch)
	{
		clear_d<Unit>(dy);
		for_each_img(src, rd, [&dy, ch, &param, rd](
			const vector<vector<PImage>> &src,
			int x, int y,
			PSize bs,
			PPoint p,
			PPoint cp
			)
		{
			if (src.size() - 1 == y)return;

			auto img1 = src[y][x];
			auto img2 = src[y + 1][x];

			img1 = clip_image(img1, p.x, img1.height - rd * 2, bs.width, rd * 2);
			img2 = clip_image(img2, p.x, 0, bs.width, rd * 2);

			auto dst = clip_image(dy, cp.x, cp.y + bs.height - rd, bs.width, rd * 2);

			y_d_p<Unit>(dst, img1, img2, ch, param.max_grandient);
		});
	},
		format,
		param
		);

}


void poisson_stiching(
	const bear::PImage &dst,
	const PStichingVectorSrc &src,
	unsigned int rd,
	unsigned int format,
	PStichingParam param)
{
	assert(
		n_channel(dst) == FI_BPP(format) &&
		8 == depth(dst));

	for_each_img(src.src, rd, [&dst](
		const vector<vector<PImage>> &src,
		int x,int y,
		PSize bs,
		PPoint p,
		PPoint cp
		) 
	{
		auto img = src[y][x];

		PImage s = clip_image(img, PRect(p, bs));

		PImage d = clip_image(dst, PRect(cp, bs));

		copy(d, s);
	});

	if (param.float_precision)
	{
		_poisson_stiching_m<float>(dst, src.src, rd, format, param);
	}
	else
	{
		_poisson_stiching_m<unsigned short>(dst, src.src, rd, format, param);
	}
}