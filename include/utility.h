#ifndef __UTILITY_HPP
#define __UTILITY_HPP

#include "../../bear/include/dynamic_image.h"
#include "../../bear/include/functor.h"

using debug_callback_t = bear::functor<void, bear::dynamic_image_ptr, bear::const_string_ptr>;

void set_debug_callback(debug_callback_t callback);

void call_debug_callback(bear::dynamic_image_ptr img, bear::const_string_ptr flag = bear::const_string_ptr());

template<int SFT, typename T>
inline T round_shift(T v)
{
	v = v + (1 << (SFT - 1));

	return v >> (SFT);
}

inline int limiteU16(int lum)
{
	if (lum & 0xffff0000)
	{
		if (lum < 0)return 0;
		return 65535;
	}
	return lum;
}


inline int limiteU8(int lum)
{
	if (lum & 0xffffff00)
	{
		if (lum < 0)return 0;
		return 255;
	}
	return lum;
}

#endif