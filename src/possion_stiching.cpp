#include "../include/possion_stiching.h"
#include "../include/possion_solver.h"

#include <cmath>

#include <algorithm>
#include <numeric>
#include <vector>

#include "../include/utility.hpp"

#include "../include/possion_stiching_dif.hpp"

#include "../../bear/include/ptr_algorism.h"

using namespace bear;
using namespace std;

template<typename Image>
image<unsigned char, 1> mask_to_byte_inner(Image src)
{
	image<unsigned char, 1> ret;
	zip_to<2>([](unsigned char &r, const_array_ptr<typename Image::elem_type> s)
	{
		r = s[0] != 0 ? 255 : 0;
	}, ret, src);
	return ret;
}

static image<unsigned char,1> mask_to_byte(const const_dynamic_image_ptr &mask)
{
	if (mask.channel_size() == 1 && mask.elm_size() == 1)
	{
		return image<unsigned char, 1>();
	}

	switch (mask.elm_size())
	{
	case 1:
		return mask_to_byte_inner(const_tensor_ptr<unsigned char, 3>(mask));
	case 2:
		return mask_to_byte_inner(const_tensor_ptr<unsigned short, 3>(mask));
	case 4:
		return mask_to_byte_inner(const_tensor_ptr<unsigned int, 3>(mask));
	case 8:
		return mask_to_byte_inner(const_tensor_ptr<unsigned long long, 3>(mask));
	default:
		throw bear_exception(exception_type::other_error, "wrong mask type!");
	}
}

template<typename Image, typename ConstImage>
void mask_merg(
	Image dst,
	ConstImage src1,
	ConstImage src2,
	const_image_ptr<unsigned char, 1> mask)
{
	using Unit = typename Image::elem_type;

	zip_to<2>([](
		array_ptr<Unit> d,
		const_array_ptr<Unit> s1,
		const_array_ptr<Unit> s2,
		unsigned char m)
	{
		if (0 == m)
		{
			copy(d, s1);
		}
		else
		{
			copy(d, s2);
		}
	}, dst, src1, src2, mask);
}

template<typename Image, typename Ds>
static void stiching(Image dst, Ds ds, unsigned int ch)
{
	using Unit = typename Ds::elm_type;

	zip_to<2>([](
		array_ptr<typename Image::elem_type> d,
		Unit s)
	{
		_ZERO_PS<Unit>::to_unit(d[ch], s);
	}, dst, ds);
}


template<typename Unit, typename Image, typename Dx, typename Dy>
static void _poisson_stiching_inner(
	Image dst,
	Dx && get_dx,
	Dy && get_dy,
	unsigned int format,
	PStichingParam param)
{
	image<Unit, 1> dx(size(dst), 1, sizeof(Unit) << 3);
	image<Unit, 1> dy(size(dst), 1, sizeof(Unit) << 3);
	image<Unit, 1> ds(size(dst), 1, sizeof(Unit) << 3);

	int chl[3];

	chl[0] = FI_RED(format);
	chl[1] = FI_GREEN(format);
	chl[2] = FI_BLUE(format);

	for (int i = 0; i < 3; ++i)
	{
		std::forward<Dx>(get_dx)(dx, chl[i]);
		std::forward<Dy>(get_dy)(dy, chl[i]);

		dxy_poisson_solver(ds, dx, dy, param.iteration_time, param.base_level);

		stiching(dst, ds, chl[i]);
	}
}

template<typename Unit, typename Image, typename ConstImage>
static void _poisson_stiching_merged(
	Image dst,
	ConstImage src,
	const_image_ptr<unsigned char, 1> _mask,
	unsigned int format,
	PStichingParam param)
{
	copy(dst, src);

	assert(src.width() == mask.width() && src.height() == mask.height());

	_poisson_stiching_inner<Unit>(dst,
		[&](image_ptr<Unit, 1> dx, int ch)
	{
			x_d<Unit>(dx, src, ch, mask, _ZERO_PS<Unit>());
	},
		[&](image_ptr<Unit, 1> dy, int ch)
	{
			y_d<Unit>(dy, src, ch, mask, _ZERO_PS<Unit>());
	},
		format,
		param
	);
}

void poisson_stiching_merged(
	dynamic_image_ptr dst,
	const_dynamic_image_ptr src,
	const_dynamic_image_ptr _mask,
	unsigned int format,
	PStichingParam param)
{
	assert(
		dst.channel_size() == FI_BPP(format) &&
		src.channel_size() == FI_BPP(format) &&
		(8 == dst.elm_size() || 16 == dst.elm_size()) &&
		(8 == src.elm_size() || 16 == src.elm_size()) &&
		src.elm_size() == dst.elm_size() &&
		src.width() == dst.width() &&
		src.height() == dst.height());


	auto maskbuf = mask_to_byte(_mask);
	const_image_ptr<unsigned char, 1> mask;
	if (maskbuf.size())
	{
		mask = maskbuf;
	}
	else
	{
		mask = const_image_ptr<unsigned char, 1>(_mask);
	}

	if (param.float_precision)
	{
		if (8 == dst.elm_size())
		{
			_poisson_stiching_merged<float>(
				tensor_ptr<unsigned char, 3>(dst),
				const_tensor_ptr<unsigned char, 3>(src),
				mask, format, param);
		}
		else
		{
			_poisson_stiching_merged<float>(
				tensor_ptr<unsigned short, 3>(dst),
				const_tensor_ptr<unsigned short, 3>(src),
				mask, format, param);

		}
	}
	else
	{
		if (8 == dst.elm_size())
		{
			_poisson_stiching_merged<unsigned short>(
				tensor_ptr<unsigned char, 3>(dst),
				const_tensor_ptr<unsigned char, 3>(src),
				mask, format, param);
		}
		else
		{
			_poisson_stiching_merged<unsigned short>(
				tensor_ptr<unsigned short, 3>(dst),
				const_tensor_ptr<unsigned short, 3>(src),
				mask, format, param);

		}
	}
}


void poisson_stiching(
	dynamic_image_ptr dst,
	const_dynamic_image_ptr src1,
	const_dynamic_image_ptr src2,
	const_dynamic_image_ptr _mask,
	unsigned int format,
	PStichingParam param)
{
	auto maskbuf = mask_to_byte(_mask);
	const_image_ptr<unsigned char, 1> mask;
	if (maskbuf.size())
	{
		mask = maskbuf;
	}
	else
	{
		mask = const_image_ptr<unsigned char, 1>(_mask);
	}

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


static void covariance(float &cv, float &max_cv, float &mse, PImage img1, PImage img2, unsigned int format)
{
	int chl[] =
	{
		FI_RED(format),
		FI_GREEN(format),
		FI_BLUE(format),
	};

	int bpp = FI_BPP(format);

	unsigned long long eab = 0;
	unsigned long long eaa = 0;
	unsigned long long ebb = 0;
	unsigned long long ea = 0;
	unsigned long long eb = 0;
	unsigned long long eapb = 0;

	for (int y = 0; y < img1.height; ++y)
	{
		unsigned char * r1 = scanline(img1, y);
		unsigned char * r2 = scanline(img2, y);

		for (int x = 0; x < img1.width; ++x)
		{
			for (int i = 0; i < 3; ++i)
			{
				unsigned int a = r1[chl[i]];
				unsigned int b = r2[chl[i]];
				eab += a * b;
				eaa += a * a;
				ebb += b * b;
				ea += a;
				eb += b;

				int apb = a - b;

				eapb += apb * apb;
			}

			r1 += bpp;
			r2 += bpp;
		}
	}

	long long cab = eab - ea * eb;
	unsigned long long caa = eaa - ea * ea;
	unsigned long long cbb = ebb - eb * eb;

	float cad2 = (float)img1.width * (float)img1.width * (float)img1.height * (float)img1.height;

	cad2 = 1.0f / cad2;

	cv = (float)cab * cad2;
	max_cv = sqrt((float)caa) * sqrt((float)cbb) * cad2;
	mse = (float)eapb * cad2;
}

inline float error_tt(PImage img1, PImage img2, unsigned int format, float th_cv, float th_mse)
{
	float cv, max_cv, mse;

	covariance(cv, max_cv, mse, img1, img2, format);

	float fcv = atan((max_cv - cv) / th_cv) * (2.0f / 3.141593653f);

	float fmse = atan(mse / th_mse) * (2.0f / 3.141593653f);

	return std::max(fcv, fmse);
}

void poisson_stiching_check(
	std::vector<bear::PPoint> error_block,
	const std::vector<bear::PPoint> eliminate,
	const PStichingVectorSrc &src,
	unsigned int rd,
	unsigned int format,
	float th_cv,
	float th_mse)
{
	th_cv = th_cv * th_cv;
	th_mse = th_mse * th_mse;

	int bw = (int)src.src[0].size();
	int bh = (int)src.src.size();

	Image dx(PSize(bw - 1, bh), 1, 32);
	Image dy(PSize(bw, bh - 1), 1, 32);

	for_each_img(src.src, rd, [=, &dx, &dy](
		const vector<vector<PImage>> &src,
		int x, int y,
		PSize bs,
		PPoint p,
		PPoint cp
		)
	{
		if (bw - 1 != x)
		{
			auto img1 = src[y][x];
			auto img2 = src[y][x + 1];

			img1 = clip_image(img1, img1.width - rd * 2, p.y, rd * 2, bs.height);
			img2 = clip_image(img2, 0, p.y, rd * 2, bs.height);

			*(float *)pick_pixel(dx, x, y) = error_tt(img1, img2, format, th_cv, th_mse);
		}

		if (bh - 1 != y)
		{
			auto img1 = src[y][x];
			auto img2 = src[y + 1][x];

			img1 = clip_image(img1, p.x, img1.height - rd * 2, bs.width, rd * 2);
			img2 = clip_image(img2, p.x, 0, bs.width, rd * 2);

			*(float *)pick_pixel(dy, x, y) = error_tt(img1, img2, format, th_cv, th_mse);
		}
	});

	enum
	{
		NORMAL_BLOCK,
		ERROR_BLOCK,
		ELIMINATE_BLOCK,
	};

	Image flag(PSize(bw, bh), 1, 32);

	fill(flag, (unsigned int)NORMAL_BLOCK);

	for (int i = 0; i < (int)eliminate.size(); ++i)
	{
		assert((unsigned int)eliminate[i].x < (unsigned int)bw && (unsigned int)eliminate[i].y < (unsigned int)bh);

		*(unsigned int *)pick_pixel(flag, eliminate[i].x, eliminate[i].y) = ELIMINATE_BLOCK;
	}

	for(;;)
	{
		int mx = 0;
		int my = 0;

		float me = 0.0f;

		for (int y = 0; y < bh; ++y)
		{
			for (int x = 0; x < bw; ++x)
			{
				if (NORMAL_BLOCK != *(unsigned int *)pick_pixel(flag, x, y))continue;

				float ce = 0;

				if (x > 0 && ERROR_BLOCK != *(unsigned int *)pick_pixel(flag, x - 1, y))
				{
					ce += *(float *)pick_pixel(dx, x - 1, y);
				}

				if (x < bw - 1 && ERROR_BLOCK != *(unsigned int *)pick_pixel(flag, x + 1, y))
				{
					ce += *(float *)pick_pixel(dx, x, y);
				}

				if (y > 0 && ERROR_BLOCK != *(unsigned int *)pick_pixel(flag, x, y - 1))
				{
					ce += *(float *)pick_pixel(dy, x, y - 1);
				}

				if (y < bh - 1 && ERROR_BLOCK != *(unsigned int *)pick_pixel(flag, x, y + 1))
				{
					ce += *(float *)pick_pixel(dy, x, y);
				}

				if (ce > me)
				{
					me = ce;
					mx = x;
					my = y;
				}
			}
		}

		if (me < 0.5f)break;

		*(unsigned int *)pick_pixel(flag, mx, my) = ERROR_BLOCK;
		error_block.push_back(PPoint(mx, my));
	}
}