#ifndef _POSSION_SOLVER_H
#define _POSSION_SOLVER_H

#include "../../bear/include/dynamic_image.h"

void dxy_poisson_solver(
	bear::dynamic_image_ptr dst,
	bear::dynamic_image_ptr dx,
	bear::dynamic_image_ptr dy,
	unsigned int iteration_time = 10,
	int base_level = 0);


void dxy_poisson_solver(
	bear::dynamic_image_ptr dst,
	bear::dynamic_image_ptr dx,
	bear::dynamic_image_ptr dy,
	bear::dynamic_image_ptr xborder,
	bear::dynamic_image_ptr yborder,
	unsigned int iteration_time = 10,
	int base_level = 0);

#endif