#ifndef _POSSION_SOLVER_H
#define _POSSION_SOLVER_H

#include "image.h"

void dxy_poisson_solver(
	const bear::PImage &dst,
	const bear::PImage &dx,
	const bear::PImage &dy,
	unsigned int iteration_time = 10,
	int base_level = 0);

#endif