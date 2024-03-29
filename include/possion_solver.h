#ifndef _POSSION_SOLVER_H
#define _POSSION_SOLVER_H

#include <bear/dynamic_image.h>
#include <functional>

extern std::function<void(bear::const_dynamic_image_ptr)> image_debug;

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

void dxy_poisson_solver(
	bear::dynamic_image_ptr dst,
	bear::dynamic_image_ptr dx,
	bear::dynamic_image_ptr dy,
	bear::dynamic_image_ptr constrain,
	unsigned int iteration_time = 10,
	int base_level = 0);

#endif