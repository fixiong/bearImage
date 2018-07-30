#ifndef _POSSION_STICHING_H
#define _POSSION_STICHING_H

#include "image.h"

void poisson_stiching(
	const bear::PImage &dst,
	const bear::PImage &src1,
	const bear::PImage &src2,
	const bear::PImage &mask,
	unsigned int format,
	int max_grandient = -1,
	unsigned int iteration_time = 10,
	int base_level = 0);

void poisson_stiching_merged(
	const bear::PImage &dst,
	const bear::PImage &src,
	const bear::PImage &mask,
	unsigned int format,
	int max_grandient = -1,
	unsigned int iteration_time = 10,
	int base_level = 0);

#endif