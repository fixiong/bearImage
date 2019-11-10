#ifndef _POSSION_STICHING_H
#define _POSSION_STICHING_H

#include <vector>
#include <algorithm>
#include "utility.h"
#include <bear/dynamic_image.h>
#include <bear/tensor.h>

enum PossionConstrain
{
	PossionNoConstrain,
	PossionEstimateConstrain,
	PossionPanoramaConstrain,
	PossionPanoramaBorderConstrain,
};

struct PStichingParam
{
	int max_grandient = -1;
	unsigned int iteration_time = 10;
	int base_level = 0;
	bool float_precision = false;

	bool edge_restriction = false;
	float edge_smooth = 1.0f;

	bear::const_dynamic_image_ptr panorama_border;

	PossionConstrain constrain = PossionNoConstrain;
};

void poisson_stiching_merged(
	bear::dynamic_image_ptr dst,
	bear::const_dynamic_image_ptr src,
	bear::const_dynamic_image_ptr mask,
	PStichingParam param);


void poisson_stiching(
	bear::dynamic_image_ptr dst,
	bear::const_dynamic_image_ptr src1,
	bear::const_dynamic_image_ptr src2,
	bear::const_dynamic_image_ptr mask,
	PStichingParam param);


struct PStichingVectorSrc
{
	template<typename Src>
	PStichingVectorSrc(Src &&_src)
	{
		int h = (int)_src.size();
		if (h <= 0)
		{
			throw bear::bear_exception(bear::exception_type::pointer_outof_range, "src is empty!");
		}
		int w = (int)_src[0].size();

		src = bear::tensor<bear::const_dynamic_image_ptr, 2>(h, w);

		for (int y = 0; y < h; ++y)
		{
			if (_src[y].size() != w)
			{
				throw bear::bear_exception(bear::exception_type::size_different, "src size inconsistence!");

			}
			for (int x = 0; x < w; ++x)
			{
				src[y][x] = std::forward<Src>(_src)[y][x];
			}
		}
	}

	bear::tensor<bear::const_dynamic_image_ptr, 2> src;
};

void poisson_stiching(
	bear::dynamic_image_ptr dst,
	const PStichingVectorSrc &src,
	size_t redundance,
	PStichingParam param = PStichingParam());


void poisson_stiching_check(
	std::vector<bear::image_point> error_block,
	bear::const_array_ptr<bear::image_point> eliminate,
	const PStichingVectorSrc &src,
	size_t rd,
	float th_cv,
	float th_mse);

void make_panorama_border(
	bear::dynamic_image_ptr border,
	bear::const_array_ptr<bear::const_dynamic_image_ptr> src);

#endif