#ifndef _POSSION_STICHING_H
#define _POSSION_STICHING_H

#include <vector>
#include <algorithm>
#include "../../bear/include/dynamic_image.h"

struct PStichingParam
{
	int max_grandient = -1;
	unsigned int iteration_time = 10;
	int base_level = 0;
	bool float_precision = false;

	bool edge_restriction = false;
	float edge_smooth = 1.0f;

};


struct MakeMask
{
	struct AutoLimite
	{
		AutoLimite(float _ignore = 0.2f, float _relax = 1.1f) :
			enable(true), ignore(_ignore), relax(_relax) {}

		bool enable;
		float ignore;
		float relax;

		static AutoLimite make_default()
		{
			AutoLimite ret;
			ret.enable = false;
			return ret;
		}
	};

	MakeMask(const bear::const_dynamic_image_ptr &_mask) :mask(_mask) {}

	template<typename Cnt>
	MakeMask(
		const Cnt &_w_border,
		const Cnt &_h_border,
		const AutoLimite &param = AutoLimite::make_default())
	{
		std::for_each(
			_w_border.begin(),
			_w_border.end(), 
			[&](typename Cnt::value_type v) 
		{
			w_border.push_back((unsigned int)v);
		});

		std::for_each(
			_h_border.begin(),
			_h_border.end(),
			[&](typename Cnt::value_type v)
		{
			h_border.push_back((unsigned int)v);
		});

		fac = param;
	}

	bear::const_dynamic_image_ptr mask;
	std::vector<unsigned int> w_border, h_border;
	AutoLimite fac;
};


void poisson_stiching(
	const bear::dynamic_image_ptr &dst,
	const bear::const_dynamic_image_ptr &src1,
	const bear::const_dynamic_image_ptr &src2,
	const bear::const_dynamic_image_ptr &mask,
	unsigned int format,
	PStichingParam param = PStichingParam());

void poisson_stiching_merged(
	const bear::dynamic_image_ptr &dst,
	const bear::const_dynamic_image_ptr &src,
	MakeMask &&mask,
	unsigned int format,
	PStichingParam param = PStichingParam());


struct PStichingVectorSrc
{
	template<typename Src>
	PStichingVectorSrc(Src &&_src)
	{
		src.resize(_src.size());
		for (int y = 0; y < (int)_src.size(); ++y)
		{
			src[y].resize(_src[0].size());
			for (int x = 0; x < (int)_src[y].size(); ++x)
			{
				src[y][x] = std::forward<Src>(_src)[y][x];
			}
		}
	}

	std::vector<std::vector<bear::const_dynamic_image_ptr>> src;
};

void poisson_stiching(
	const bear::dynamic_image_ptr &dst,
	const PStichingVectorSrc &src,
	unsigned int redundance,
	unsigned int format,
	PStichingParam param = PStichingParam());


void poisson_stiching_check(
	std::vector<bear::image_point> error_block,
	const std::vector<bear::image_point> eliminate,
	const PStichingVectorSrc &src,
	unsigned int rd,
	unsigned int format,
	float th_cv,
	float th_mse);

#endif