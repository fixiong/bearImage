#ifndef _POSSION_SOLVER_H
#define _POSSION_SOLVER_H

#include "../../bear/include/dynamic_image.h"

void dxy_poisson_solver(
	const bear::dynamic_image_ptr &dst,
	const bear::dynamic_image_ptr &dx,
	const bear::dynamic_image_ptr &dy,
	unsigned int iteration_time = 10,
	int base_level = 0);

#endif