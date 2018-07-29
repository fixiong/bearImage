#ifndef _TYPES_H
#define _TYPES_H

#define F_RGBA (0 | (1 << 8) | (2 << 16) | (3 << 24) | (4 << 28))
#define F_BGRA (2 | (1 << 8) | (0 << 16) | (3 << 24) | (4 << 28))
#define F_ARGB (1 | (2 << 8) | (3 << 16) | (0 << 24) | (4 << 28))
#define F_ABGR (3 | (2 << 8) | (1 << 16) | (0 << 24) | (4 << 28))
#define F_RGB (0 | (1 << 8) | (2 << 16) | (15 << 24) | (3 << 28))
#define F_BGR (2 | (1 << 8) | (0 << 16) | (15 << 24) | (3 << 28))

#define FI_RED(n)   (n & 255)
#define FI_GREEN(n) ((n >> 8) & 255)
#define FI_BLUE(n)  ((n >> 16) & 255)
#define FI_ALPHA(n) ((n >> 24) & 15)
#define FI_BPP(n)   ((n >> 28) & 15)


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