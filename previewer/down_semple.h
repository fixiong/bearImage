#pragma once

#include <bear/image.h>

void down_semple(
	bear::image_ptr<unsigned char, 3> dst,
	bear::const_image_ptr<unsigned char, 3> src,
	float x_first_Pos,
	float x_fac,
	float y_first_Pos,
	float y_fac);