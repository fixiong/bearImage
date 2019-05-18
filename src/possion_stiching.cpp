#include "../include/possion_stiching.h"
#include "../include/possion_solver.h"

#include <cmath>

#include <algorithm>
#include <numeric>
#include <vector>

#include "../include/utility.h"
#include <bear/ptr_algorism.h>
#include <bear/ptr_numeric.h>

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
		//typename Image::elm_type u = 128;
		//_ZERO_PS<Unit>::to_unit(u, s);
		//d[ch] = 255 - u;

		_ZERO_PS<Unit>::to_unit(d[ch], s);

	}, dst, ds);
}

template<typename Unit, typename Array>
static void estimate_border(Array &db, bool ring)
{
	using CT = decltype(db[0] + db[0]);
	vector<CT> acc(db.size());

	acc[0] = db[0];

	for (int i = 1; i < (int)db.size(); ++i)
	{
		acc[i] = _ZERO_PS<Unit>::acc(acc[i - 1], db[i]);
	}

	if (ring)
	{
		auto er = _ZERO_PS<Unit>::run() - acc.back();

		for (int i = 0; i < (int)db.size(); ++i)
		{
			acc[i] += er * (i + 1) / (int)db.size();
		}
	}

	auto avg = accumulate(acc.begin(), acc.end(), (CT)0) / (int)db.size();

	map_function([avg](Unit &d, CT s)
	{
		d = _ZERO_PS<Unit>::limite(s - avg);
	}, db, acc);
}


template<typename Unit, typename Image, typename Dx, typename Dy>
static void _poisson_stiching_inner(
	Image dst,
	Dx && get_dx,
	Dy && get_dy,
	PStichingParam param)
{
	image<Unit, 1> dx(width(dst), height(dst));
	image<Unit, 1> dy(width(dst), height(dst));
	image<Unit, 1> ds(width(dst), height(dst));

	for (int i = 0; i < dst[0][0].size(); ++i)
	{
		std::forward<Dx>(get_dx)(dx, i);
		std::forward<Dy>(get_dy)(dy, i);

		if (param.constrain == PossionNoConstrain)
		{
			dxy_poisson_solver(ds, dx, dy, param.iteration_time, param.base_level);

		}
		else
		{
			int w = (int)width(dst);
			int h = (int)height(dst);
			image<Unit, 1> bx(2, h);
			image<Unit, 1> by(w, 2);


			if (param.constrain == PossionPanoramaConstrain)
			{
				vector<Unit> bd1(h);
				vector<Unit> bd2(h);

				for (int y = 0; y < h; ++y)
				{
					bd1[y] = dy[y][0];
					bd2[y] = dy[y][w - 1];
				}

				map_function([](Unit &v1, Unit v2) {
					v1 = (v1 + v2) / 2;
				}, bd1, bd2);

				estimate_border<Unit>(bd1, false);

				for (int y = 0; y < h; ++y)
				{
					bx[y][0] = bx[y][1] = bd1[y];
				}

				by[0].fill(bd1[0]);
				by[1].fill(bd1.back());
			}
			else
			{
				vector<Unit> bd(w * 2 + h * 2 - 4);

				copy(to_ptr(bd).clip(0, w), dx[0]);

				for (int y = 1; y < h; ++y)
				{
					bd[w - 1 + y] = dy[y][w - 1];
				}

				for (int x = 1; x < w; ++x)
				{
					bd[w + h - 2 + x] = _ZERO_PS<Unit>::minus(dx[h - 1][w - x]);
				}

				for (int y = 1; y < h - 1; ++y)
				{
					bd[w * 2 + h - 3 + y] = _ZERO_PS<Unit>::minus(dx[h - y][0]);
				}

				bd[0] = _ZERO_PS<Unit>::minus(dx[1][0]);

				estimate_border<Unit>(bd, true);

				copy(by[0], to_ptr(bd).clip(0, w));

				for (int y = 0; y < h; ++y)
				{
					bx[y][1] = bd[w - 1 + y];
				}

				for (int x = 0; x < w; ++x)
				{
					by[1][w - 1 - x] = bd[w + h - 2 + x];
				}

				for (int y = 0; y < h - 1; ++y)
				{
					bx[h - 1 - y][0] = bd[w * 2 + h - 3 + y];
				}

				bx[0][0] = bd[0];
			}
			dxy_poisson_solver(ds, dx, dy, bx, by, param.iteration_time, param.base_level);
		}


		stiching(dst, ds, i);
	}
}

template<typename Unit, typename Image, typename ConstImage>
static void _poisson_stiching_merged(
	Image dst,
	ConstImage src,
	const_image_ptr<unsigned char, 1> mask,
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
		param
	);
}

void poisson_stiching_merged(
	dynamic_image_ptr dst,
	const_dynamic_image_ptr src,
	const_dynamic_image_ptr _mask,
	PStichingParam param)
{
	assert(
		(1 == dst.elm_size() || 2 == dst.elm_size()) &&
		(1 == src.elm_size() || 2 == src.elm_size()) &&
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
		if (1 == dst.elm_size())
		{
			_poisson_stiching_merged<float>(
				tensor_ptr<unsigned char, 3>(dst),
				const_tensor_ptr<unsigned char, 3>(src),
				mask, param);
		}
		else
		{
			_poisson_stiching_merged<float>(
				tensor_ptr<unsigned short, 3>(dst),
				const_tensor_ptr<unsigned short, 3>(src),
				mask, param);

		}
	}
	else
	{
		if (1 == dst.elm_size())
		{
			_poisson_stiching_merged<unsigned short>(
				tensor_ptr<unsigned char, 3>(dst),
				const_tensor_ptr<unsigned char, 3>(src),
				mask, param);
		}
		else
		{
			_poisson_stiching_merged<unsigned short>(
				tensor_ptr<unsigned short, 3>(dst),
				const_tensor_ptr<unsigned short, 3>(src),
				mask,  param);

		}
	}
}


void poisson_stiching(
	dynamic_image_ptr dst,
	const_dynamic_image_ptr src1,
	const_dynamic_image_ptr src2,
	const_dynamic_image_ptr _mask,
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

	if (dst.elm_size() == 1)
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

	poisson_stiching_merged(dst, src1, mask, param);
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
			if (width(src) - 1 == x)return;

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
			if (height(src) - 1 == y)return;

			auto sz = image_size(bs.width, rd * 2);

			auto img1 = src[y][x];
			auto img2 = src[y + 1][x];

			img1 = clip_image(img1, image_rectangle(p.x, height(img1) - rd * 2, sz.width, sz.height));
			img2 = clip_image(img2, image_rectangle(p.x, 0, sz.width, sz.height));

			auto dst = dy.clip(image_rectangle(cp.x, cp.y + bs.height - rd, sz.width, sz.height));

			y_d_p(dst, img1, img2, ch, _ZERO_PS<Unit>());
		});
	},
		param);

}

template<typename Dst, typename Src>
static void _poisson_stiching_a(
	Dst dst,
	Src src,
	size_t rd,
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
		_poisson_stiching_m<float>(dst, src, rd, param);
	}
	else
	{
		_poisson_stiching_m<unsigned short>(dst, src, rd, param);
	}

}

static image_size check_src(const_tensor_ptr<const_dynamic_image_ptr, 2> src, size_t redundance)
{
	if (0 == to_ptr(src).total_size())
	{
		throw bear_exception(exception_type::pointer_outof_range, "empty src!");
	}

	auto es = src[0][0].elm_size();

	if (es != 1 && es != 2)
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
			if (0 == y)
			{
				ws[x] = img.width();
			}

			if (ws[x] != img.width() ||
				hs[y] != img.height() ||
				es != img.elm_size())
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

	ret.height -= redundance * (h - 1) * 2;
	ret.width -= redundance * (w - 1) * 2;

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
	PStichingParam param)
{
	if (dst.size() != check_src(src.src, rd) || dst.elm_size() != src.src[0][0].elm_size())
	{
		throw bear_exception(exception_type::size_different, "src dst size different!");
	}

	if (1 == dst.elm_size())
	{
		_poisson_stiching_a(
			tensor_ptr<unsigned char, 3>(dst),
			to_ptr(map_function([] (const const_dynamic_image_ptr & img) 
				-> wrapper<const_tensor_ptr<unsigned char, 3>> 
		{
			return const_tensor_ptr<unsigned char, 3>(img);
		}, src.src)),
			rd,param);
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
			rd, param);
	}
}

template<typename Image>
static void covariance(
	float &cv,
	float &max_cv,
	float &mse,
	Image img1,
	Image img2)
{
	using Unit = typename Image::elm_type;

	unsigned long long eab = 0;
	unsigned long long eaa = 0;
	unsigned long long ebb = 0;
	unsigned long long ea = 0;
	unsigned long long eb = 0;
	unsigned long long eapb = 0;

	zip_to<3>([&](Unit r1, Unit r2 )
	{
		unsigned int a = r1;
		unsigned int b = r2;
		eab += a * b;
		eaa += a * a;
		ebb += b * b;
		ea += a;
		eb += b;

		int apb = a - b;

		eapb += apb * apb;
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
	float th_cv,
	float th_mse)
{
	float cv, max_cv, mse;

	covariance(cv, max_cv, mse, img1, img2);

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
	float th_cv,
	float th_mse)
{
	th_cv = th_cv * th_cv;
	th_mse = th_mse * th_mse;

	int bw = (int)width(src);
	int bh = (int)height(src);

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

			dx[y][x] = error_tt(img1, img2, th_cv, th_mse);
		}

		if (bh - 1 != y)
		{
			auto img1 = src[y][x];
			auto img2 = src[y + 1][x];

			img1 = clip_image(img1, image_rectangle(p.x, height(img1) - rd * 2, bs.width, rd * 2));
			img2 = clip_image(img2, image_rectangle(p.x, 0, bs.width, rd * 2));

			dy[y][x] = error_tt(img1, img2, th_cv, th_mse);
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
	float th_cv,
	float th_mse)
{
	check_src(src.src, rd);

	if (1 == src.src[0][0].elm_size())
	{
		_poisson_stiching_check(
			error_block,
			eliminate,
			to_ptr(map_function([](const const_dynamic_image_ptr & img)
				-> wrapper<const_tensor_ptr<unsigned char, 3>>
		{
			return const_tensor_ptr<unsigned char, 3>(img);
		}, src.src)),
			rd, th_cv, th_mse);
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
			rd, th_cv, th_mse);
	}
}