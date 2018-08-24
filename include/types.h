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
		PSize() = default;
		PSize(unsigned int _width, unsigned int _height) :
			width(_width), height(_height) {}

		int width = 0;
		int height = 0;

		bool operator == (const PSize &oth) const
		{
			return width == oth.width && height == oth.height;
		}
	};


	struct PPoint
	{
		PPoint() = default;
		PPoint(int _x,int _y) :
			x(_x), y(_y) {}
		int x = 0;
		int y = 0;

		bool operator == (const PPoint &oth) const
		{
			return x == oth.x && y == oth.y;
		}
	};

	struct PRect
	{
		PRect() = default;
		PRect(int _x,int _y,unsigned int _width, unsigned int _height) :
			x(_x),y(_y),
			width(_width), height(_height) {}


		PRect(PPoint pos,PSize size) :
			x(pos.x), y(pos.y),
			width(size.width), height(size.height) {}

		int x = 0, y = 0, width = 0, height = 0;

		PSize size() const
		{
			return PSize(width, height);
		}

		PPoint pos() const
		{
			return PPoint(x, y);
		}

		bool operator == (const PRect &oth) const
		{
			return size() == oth.size() && pos() == oth.pos();
		}
	};
}


#endif