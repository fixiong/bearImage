#include "image.h"

#include <assert.h>
#include <memory.h>

namespace bear
{
	void copy(const PImage &dst, const PImage &src)
	{
		unsigned int cl = src.width * src.depth * src.n_channel / 8;
		assert(dst.height >= src.height &&
			dst.width * dst.depth * dst.n_channel / 8 >= cl);

		if (src.data == dst.data)return;

		for (int i = 0; i<src.height; ++i)
		{
			memmove(scanline(dst, i), scanline(src, i), cl);
		}

	}

	void zero(const PImage &img)
	{
		for (int y = 0; y < height(img); ++y)
		{
			memset(scanline(img, y), 0, img.width * (img.depth >> 3) * img.n_channel);
		}
	}

	void zero(const PImage &img,unsigned int ch)
	{
		assert(ch < img.n_channel);

		for (int y = 0; y < height(img); ++y)
		{
			unsigned b = (depth(img) >> 3);

			unsigned int bpp = b * n_channel(img);

			unsigned char * row = scanline(img, y) + ch * b;
			for (int x = 0; x < width(img); ++x)
			{
				for (int i = 0; i < b; ++i)
				{
					row[i] = 0;
				}

				row += bpp;
			}
		}
	}


	void assert_range(const PImage &img, void * ptr)
	{
		assert(
			(unsigned char *)ptr >= scanline(img, 0) ||
			(unsigned char *)ptr < scanline(img, img.height)
		);
	}

	PImage clip_image(const PImage &img, int x_offset, int y_offset, int width, int height)
	{
		assert(
			x_offset + width <= img.width ||
			y_offset + height <= img.height
		);

		PImage ret;
		ret.data = pick_pixel(img, x_offset, y_offset);
		ret.depth = img.depth;
		ret.n_channel = img.n_channel;
		ret.width = width;
		ret.height = height;
		ret.width_step = img.width_step;
		return ret;
	}


	Image::Image(Image &&other)
	{
		width = other.width;
		height = other.height;
		depth = other.depth;
		n_channel = other.n_channel;
		width_step = other.width_step;
		data = other.data;
		buf = other.buf;

		other.buf = 0;
	}

	Image & Image::operator =(Image&& other)
	{
		destruct();

		width = other.width;
		height = other.height;
		depth = other.depth;
		n_channel = other.n_channel;
		width_step = other.width_step;
		data = other.data;
		buf = other.buf;

		other.buf = 0;

		return *this;
	}

	template<typename T>
	inline T align_m(T p, unsigned int align)
	{
		if (align)
		{
			unsigned int m = p % align;

			if (m)
			{
				return p + align - m;
			}
		}
		return p;
	}

	void Image::construct(
		unsigned int _width,
		unsigned int _height,
		unsigned int _n_channel,
		unsigned int _depth,
		unsigned int align)
	{
		destruct();

		static unsigned int min_align = 16;
		width = _width;
		height = _height;
		n_channel = _n_channel;
		depth = _depth;
		width_step = align_m(width * n_channel * depth / 8, align);

		if (align <= min_align)
		{
			data = buf = new unsigned char[width_step * height];

			if (!((CompatibleType)buf % align))return;

			min_align = align >> 1;

			delete[] buf;
		}

		buf = new unsigned char[width_step * height + align];

		data = (unsigned char *)align_m((CompatibleType)buf, align);
	}

	void Image::destruct()
	{
		if (buf)
		{
			delete[] buf;
			buf = 0;
		}
	}
}