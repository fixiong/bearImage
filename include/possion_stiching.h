#ifndef _POSSION_STICHING_H
#define _POSSION_STICHING_H

#include "image.h"

struct PStichingParam
{
	PStichingParam()
	{
		max_grandient = -1;
		iteration_time = 10;
		base_level = 0;
		float_precision = false;
	}

	int max_grandient;
	unsigned int iteration_time;
	int base_level;
	bool float_precision;
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
	const bear::PImage &mask,
	unsigned int format,
	PStichingParam param = PStichingParam());

#endif