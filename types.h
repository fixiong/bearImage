#ifndef _TYPES_H
#define _TYPES_H

namespace bear
{
	template<unsigned int v>
	struct __COMPATIBLETYPE
	{
		typedef unsigned int Type;
	};

	template<>
	struct __COMPATIBLETYPE<8>
	{
		typedef unsigned long long Type;
	};

	typedef typename __COMPATIBLETYPE<sizeof(void *)>::Type CompatibleType;

	struct _PImage
	{
		int width;
		int height;
		unsigned int depth;
		unsigned int n_channel;
		unsigned int width_step;
		unsigned char * data;

	};

	struct PSize
	{
		PSize(unsigned int _width, unsigned int _height) :
			width(_width), height(_height) {}
		unsigned int width;
		unsigned int height;
	};
}


#endif