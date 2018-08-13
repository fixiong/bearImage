#ifndef _POSSION_STICHING_H
#define _POSSION_STICHING_H

#include "image.h"
#include <vector>
#include <algorithm>

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

	MakeMask(const bear::PImage &_mask) :mask(_mask) {}

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

	bear::PImage mask;
	std::vector<unsigned int> w_border, h_border;
	AutoLimite fac;
};


void poisson_stiching(
	const bear::PImage &dst,
	const bear::PImage &src1,
	const bear::PImage &src2,
	const bear::PImage &mask,
	unsigned int format,
	PStichingParam param = PStichingParam());

void poisson_stiching_merged(
	const bear::PImage &dst,
	const bear::PImage &src,
	MakeMask &&mask,
	unsigned int format,
	PStichingParam param = PStichingParam());

#endif