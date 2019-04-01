#include "../include/utility.h"

static debug_callback_t callback;

void set_debug_callback(debug_callback_t _callback)
{
	callback = _callback;
}

void call_debug_callback(bear::dynamic_image_ptr img, bear::const_string_ptr flag)
{
	auto cb = callback;
	if (cb)
	{
		cb(img, flag);
	}
}
