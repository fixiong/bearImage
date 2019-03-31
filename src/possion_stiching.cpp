#include "../include/possion_stiching.h"
#include "../include/possion_solver.h"

#include <cmath>

#include <algorithm>
#include <numeric>
#include <vector>

#include "../include/utility.h"
#include "../../bear/include/ptr_algorism.h"

#include "possion_stiching_dif.hpp"


using namespace bear;
using namespace std;

template<typename Image>
image<unsigned char, 1> mask_to_byte_inner(Image src)
{
	image<unsigned char, 1> ret;
	zip_to<2>([](unsigned char &r, const_array_ptr<typename Image::elm_type> s)
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
	using Unit = typename Image::elm_type;

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

	zip_to<2>([ch](
		array_ptr<typename Image::elm_type> d,
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
	image<Unit, 1> dx(width(dst), height(dst));
	image<Unit, 1> dy(width(dst), height(dst));
	image<Unit, 1> ds(width(dst), height(dst));

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
	const_image_ptr<unsigned char, 1> mask,
	unsigned int format,
	PStichingParam param)
{
	copy(dst, src);

	assert(width(src) == mask.width() && height(src) == mask.height());

	_poisson_stiching_inner<Unit>(dst,
		[&](image_ptr<Unit, 1> dx, int ch)
	{
			x_d(dx, src, ch, mask, _ZERO_PS<Unit>());
	},
		[&](image_ptr<Unit, 1> dy, int ch)
	{
			y_d(dy, src, ch, mask, _ZERO_PS<Unit>());
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

	if (dst.elm_size() == 8)
	{
		mask_merg(
			tensor_ptr<unsigned char, 3>(dst),
			const_tensor_ptr<unsigned char, 3>(src1),
			const_tensor_ptr<unsigned char, 3>(src2),
			mask);
	}
	else
	{
		mask_merg(
			tensor_ptr<unsigned short, 3>(dst),
			const_tensor_ptr<unsigned short, 3>(src1),
			const_tensor_ptr<unsigned short, 3>(src2),
			mask);
	}

	poisson_stiching_merged(dst, src1, mask, format, param);
}


template<typename Src, typename C>
static void for_each_img(Src src,size_t rd, C && c)
{
	auto nw = width(src);
	auto nh = height(src);

	size_t cpy = 0;
	for (size_t y = 0; y < nh; ++y)
	{

		auto rh = height(src[y][0]);
		auto bh = rh;
		size_t py = 0;

		if (0 != y)
		{
			bh -= rd;
			py += rd;
		}

		if (nh - 1 != y) bh -= rd;

		size_t cpx = 0;
		for (size_t x = 0; x < nw; ++x)
		{

			auto rw = width(src[y][x]);
			auto bw = rw;
			size_t px = 0;

			if (0 != x)
			{
				bw -= rd;
				px += rd;
			}

			if (nw - 1 != x) bw -= rd;

			forward<C>(c)(
				x, y,
				image_size(bw, bh),
				image_point(px, py),
				image_point(cpx, cpy));

			cpx += bw;
		}

		cpy += bh;
	}
}

template<typename Unit, typename Dst, typename Src>
static void _poisson_stiching_m(
	Dst dst,
	Src src,
	size_t rd,
	unsigned int format,
	PStichingParam param)
{

	_poisson_stiching_inner<Unit>(
		dst,
		[&src, &param, rd](image_ptr<Unit,1> dx, int ch)
	{
		dx.fill(_ZERO_PS<Unit>::run());
		for_each_img(src, rd, [&src, &dx, ch, &param, rd](
			size_t x, size_t y,
			image_size bs,
			image_point p,
			image_point cp)
		{
			if (src[y].size() - 1 == x)return;

			auto sz = image_size(rd * 2, bs.height);

			auto img1 = src[y][x];
			auto img2 = src[y][x + 1];

			img1 = clip_image(img1, image_rectangle(width(img1) - rd * 2, p.y, sz.width, sz.height));
			img2 = clip_image(img2, image_rectangle(0, p.y, sz.width, sz.height));

			auto dst = dx.clip(image_rectangle(cp.x + bs.width - rd, cp.y, sz.width, sz.height));

			x_d_p(dst, img1, img2, ch, _ZERO_PS<Unit>());

		});
	},
		[&src, &param, rd](image_ptr<Unit, 1> dy, int ch)
	{
		dy.fill(_ZERO_PS<Unit>::run());;
		for_each_img(src, rd, [&src, &dy, ch, &param, rd](
			size_t x, size_t y,
			image_size bs,
			image_point p,
			image_point cp)
		{
			if (src.size() - 1 == y)return;

			auto sz = image_size(bs.width, rd * 2);

			auto img1 = src[y][x];
			auto img2 = src[y + 1][x];

			img1 = clip_image(img1, image_rectangle(p.x, height(img1) - rd * 2, sz.width, sz.height));
			img2 = clip_image(img2, image_rectangle(p.x, 0, sz.width, sz.height));

			auto dst = dy.clip(image_rectangle(cp.x, cp.y + bs.height - rd, sz.width, sz.height));

			y_d_p(dst, img1, img2, ch, _ZERO_PS<Unit>());
		});
	},
		format,
		param);

}

template<typename Dst, typename Src>
static void _poisson_stiching_a(
	Dst dst,
	Src src,
	size_t rd,
	unsigned int format,
	PStichingParam param)
{
	for_each_img(src, rd, [&dst, &src](
		size_t x, size_t y,
		image_size bs,
		image_point p,
		image_point cp
		)
	{
		auto img = src[y][x];

		auto s = clip_image(img, image_rectangle{ p, bs });
		auto d = clip_image(dst, image_rectangle{ cp, bs });

		copy(d, s);
	});

	if (param.float_precision)
	{
		_poisson_stiching_m<float>(dst, src, rd, format, param);
	}
	else
	{
		_poisson_stiching_m<unsigned short>(dst, src, rd, format, param);
	}

}

static image_size check_src(const_tensor_ptr<const_dynamic_image_ptr, 2> src, size_t redundance, unsigned int format)
{
	if (0 == to_ptr(src).total_size())
	{
		throw bear_exception(exception_type::pointer_outof_range, "empty src!");
	}

	auto es = src[0][0].elm_size();

	if (es != 8 && es != 16)
	{
		throw bear_exception(exception_type::other_error, "only support to 24bit or 48bit image!");
	}

	auto h = height(src);
	auto w = width(src);

	vector<size_t> hs(h);
	vector<size_t> ws(w);

	for (size_t y = 0; y < h; ++y)
	{
		hs[y] = src[y][0].height();
		for (size_t x = 0; x < w; ++x)
		{
			auto img =src[y][x];
			if (0 == x)
			{
				ws[x] = img.width();
			}

			if (ws[x] != img.width() ||
				hs[y] != img.height() ||
				es != img.elm_size() ||
				FI_BPP(format) != img.channel_size())
			{
				throw bear_exception(exception_type::size_different, "src size inconsistence!");
			}
		}
	}

	image_size ret(0,0);

	to_ptr(hs).for_each([&ret](size_t s) {
		ret.height += s;
	});

	to_ptr(ws).for_each([&ret](size_t s) {
		ret.width += s;
	});

	ret.height -= redundance * (h - 1);
	ret.width -= redundance * (w - 1);

	if (ret.height > 1000000 || ret.width > 1000000)
	{
		throw bear_exception(exception_type::pointer_outof_range, "wrong dst size!");
	}

	return ret;
}



void poisson_stiching(
	dynamic_image_ptr dst,
	const PStichingVectorSrc &src,
	size_t rd,
	unsigned int format,
	PStichingParam param)
{
	if (dst.size() != check_src(src.src, rd, format) || dst.elm_size() != src.src[0][0].elm_size())
	{
		throw bear_exception(exception_type::size_different, "src dst size different!");
	}

	if (8 == dst.elm_size())
	{
		_poisson_stiching_a(
			tensor_ptr<unsigned char, 3>(dst),
			to_ptr(map_function([] (const const_dynamic_image_ptr & img) 
				-> wrapper<const_tensor_ptr<unsigned char, 3>> 
		{
			return const_tensor_ptr<unsigned char, 3>(img);
		}, src.src)),
			rd,format,param);
	}
	else
	{
		_poisson_stiching_a(
			tensor_ptr<unsigned short, 3>(dst),
			to_ptr(map_function([](const const_dynamic_image_ptr & img)
				-> wrapper<const_tensor_ptr<unsigned short, 3>>
		{
			return const_tensor_ptr<unsigned short, 3>(img);
		}, src.src)),
			rd, format, param);
	}
}

template<typename Image>
static void covariance(
	float &cv,
	float &max_cv,
	float &mse,
	Image img1,
	Image img2,
	unsigned int format)
{
	using Unit = typename Image::elm_type;

	int chl[] =
	{
		FI_RED(format),
		FI_GREEN(format),
		FI_BLUE(format),
	};

	unsigned long long eab = 0;
	unsigned long long eaa = 0;
	unsigned long long ebb = 0;
	unsigned long long ea = 0;
	unsigned long long eb = 0;
	unsigned long long eapb = 0;

	zip_to<2>([&](array_ptr<Unit> r1, array_ptr<Unit> r2 )
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
	}, img1, img2);

	long long cab = eab - ea * eb;
	unsigned long long caa = eaa - ea * ea;
	unsigned long long cbb = ebb - eb * eb;

	float cad2 = (float)width(img1) * (float)width(img1) * (float)height(img1) * (float)height(img1);

	cad2 = 1.0f / cad2;

	cv = (float)cab * cad2;
	max_cv = sqrt((float)caa) * sqrt((float)cbb) * cad2;
	mse = (float)eapb * cad2;
}

template<typename Image>
inline float error_tt(
	Image img1,
	Image img2,
	unsigned int format,
	float th_cv,
	float th_mse)
{
	float cv, max_cv, mse;

	covariance(cv, max_cv, mse, img1, img2, format);

	float fcv = atan((max_cv - cv) / th_cv) * (2.0f / 3.141593653f);

	float fmse = atan(mse / th_mse) * (2.0f / 3.141593653f);

	return std::max(fcv, fmse);
}

template<typename Src>
void _poisson_stiching_check(
	vector<image_point> &error_block,
	const_array_ptr<image_point> eliminate,
	Src src,
	size_t rd,
	unsigned int format,
	float th_cv,
	float th_mse)
{
	th_cv = th_cv * th_cv;
	th_mse = th_mse * th_mse;

	int bw = (int)src[0].size();
	int bh = (int)src.size();

	image<float, 1> dx(bw - 1, bh);
	image<float, 1> dy(bw, bh - 1);

	for_each_img(src, rd, [=, &dx, &dy](
		size_t x, size_t y,
		image_size bs,
		image_point p,
		image_point cp
		)
	{
		if (bw - 1 != x)
		{
			auto img1 = src[y][x];
			auto img2 = src[y][x + 1];

			img1 = clip_image(img1, image_rectangle(width(img1) - rd * 2, p.y, rd * 2, bs.height));
			img2 = clip_image(img2, image_rectangle(0, p.y, rd * 2, bs.height));

			dx[y][x] = error_tt(img1, img2, format, th_cv, th_mse);
		}

		if (bh - 1 != y)
		{
			auto img1 = src[y][x];
			auto img2 = src[y + 1][x];

			img1 = clip_image(img1, image_rectangle(p.x, height(img1) - rd * 2, bs.width, rd * 2));
			img2 = clip_image(img2, image_rectangle(p.x, 0, bs.width, rd * 2));

			dy[y][x] = error_tt(img1, img2, format, th_cv, th_mse);
		}
	});

	enum
	{
		NORMAL_BLOCK,
		ERROR_BLOCK,
		ELIMINATE_BLOCK,
	};

	image<unsigned int, 1> flag(bw, bh);

	to_ptr(flag).fill((unsigned int)NORMAL_BLOCK);

	to_ptr(eliminate).for_each([bw,bh,&flag](image_point p)
	{
		if ((unsigned int)p.x >= (unsigned int)bw || (unsigned int)p.y >= (unsigned int)bh)
		{
			throw bear_exception(exception_type::pointer_outof_range, "point outside image!");
		}

		flag[p.y][p.x] = ELIMINATE_BLOCK;
	});

	for (;;)
	{
		int mx = 0;
		int my = 0;

		float me = 0.0f;

		for (int y = 0; y < bh; ++y)
		{
			for (int x = 0; x < bw; ++x)
			{
				if (NORMAL_BLOCK != flag[y][x])continue;

				float ce = 0;

				if (x > 0 && ERROR_BLOCK != flag[y][x - 1])
				{
					ce += dx[y][x - 1];
				}

				if (x < bw - 1 && ERROR_BLOCK != flag[y][x + 1])
				{
					ce += dx[y][x];
				}

				if (y > 0 && ERROR_BLOCK != flag[y - 1][x])
				{
					ce += dy[y - 1][x];
				}

				if (y < bh - 1 && ERROR_BLOCK != flag[y + 1][x])
				{
					ce += dy[y][x];
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

		flag[my][mx] = ERROR_BLOCK;
		error_block.push_back(image_point(mx, my));
	}
}

void poisson_stiching_check(
	vector<image_point> error_block,
	const_array_ptr<image_point> eliminate,
	const PStichingVectorSrc &src,
	size_t rd,
	unsigned int format,
	float th_cv,
	float th_mse)
{
	check_src(src.src, rd, format);

	if (8 == src.src[0][0].elm_size())
	{
		_poisson_stiching_check(
			error_block,
			eliminate,
			to_ptr(map_function([](const const_dynamic_image_ptr & img)
				-> wrapper<const_tensor_ptr<unsigned char, 3>>
		{
			return const_tensor_ptr<unsigned char, 3>(img);
		}, src.src)),
			rd, format, th_cv, th_mse);
	}
	else
	{
		_poisson_stiching_check(
			error_block,
			eliminate,
			to_ptr(map_function([](const const_dynamic_image_ptr & img)
				-> wrapper<const_tensor_ptr<unsigned short, 3>>
		{
			return const_tensor_ptr<unsigned short, 3>(img);
		}, src.src)),
			rd, format, th_cv, th_mse);
	}
}