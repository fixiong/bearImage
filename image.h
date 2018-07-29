#ifndef _IMAGE_H
#define _IMAGE_H

#include <assert.h>
#include <utility>
#include "types.h"

namespace bear
{

#ifdef CV_MAJOR_VERSION

	inline _PImage to_p_image(const IplImage &img)
	{
		_PImage ret;

		ret.width = img.width;
		ret.height = img.height;
		ret.n_channel = img.nChannels;
		ret.depth = img.depth;
		ret.width_step = img.widthStep;
		ret.data = (unsigned char *)img.imageData;

		return ret;
	}

	inline _PImage to_p_image(const IplImage * img)
	{
		return to_p_image(*img);
	}

	inline _PImage to_p_image(const CvMat &img)
	{
		_PImage ret;

		ret.width = img.width;
		ret.height = img.height;
		ret.n_channel = CV_MAT_CN(img.type);

		switch (CV_MAT_DEPTH(img.type))
		{
		case CV_8U:
		case CV_8S:
			ret.depth = 8;
			break;
		case CV_16U:
		case CV_16S:
			ret.depth = 16;
			break;
		case CV_32S:
		case CV_32F:
			ret.depth = 32;
			break;
		default:
			ret.depth = 8;
		}

		ret.width_step = img.step;
		ret.data = (unsigned char *)img.data.ptr;

		return ret;
	}

	inline _PImage to_p_image(const cv::Mat &img)
	{
		return to_p_image(CvMat(img));
	}

#endif

#ifdef QT_VERSION

	inline _PImage to_p_image(QImage &img)
	{
		_PImage ret;

		ret.width = img.width();
		ret.height = img.height();

		switch (img.format())
		{
		case QImage::Format_Mono:
		case QImage::Format_MonoLSB:
		case QImage::Format_Indexed8:
		case QImage::Format_Alpha8:
		case QImage::Format_Grayscale8:
			ret.n_channel = 1;
			ret.depth = 8;
			break;
		case QImage::Format_ARGB32:
		case QImage::Format_ARGB32_Premultiplied:
		case QImage::Format_RGBX8888:
		case QImage::Format_RGBA8888:
		case QImage::Format_RGBA8888_Premultiplied:
			ret.n_channel = 4;
			ret.depth = 8;
			break;
		case QImage::Format_RGB32:
		case QImage::Format_BGR30:
		case QImage::Format_A2BGR30_Premultiplied:
		case QImage::Format_RGB30:
		case QImage::Format_A2RGB30_Premultiplied:
			ret.n_channel = 1;
			ret.depth = 32;
			break;
		case QImage::Format_RGB16:
		case QImage::Format_RGB555:
		case QImage::Format_RGB444:
		case QImage::Format_ARGB4444_Premultiplied:
			ret.n_channel = 1;
			ret.depth = 16;
			break;
		case QImage::Format_ARGB8565_Premultiplied:
		case QImage::Format_RGB666:
		case QImage::Format_ARGB6666_Premultiplied:
		case QImage::Format_ARGB8555_Premultiplied:
		case QImage::Format_RGB888:
			ret.n_channel = 3;
			ret.depth = 8;
			break;
		default:
			ret.n_channel = 4;
			ret.depth = 8;
		}

		ret.width_step = img.bytesPerLine();
		ret.data = img.bits();

		return ret;
	}


#endif // QT_VERSION



	inline _PImage to_p_image(const _PImage &img)
	{
		return img;
	}


	struct PImage : public _PImage
	{
		template<typename T>
		PImage(T && other): _PImage(to_p_image(std::forward<T>(other))){}

		PImage(void *_data, int _width, int _height, int _width_step, int _n_channel, int _depth)
		{
			width = _width;
			height = _height;
			depth = _width_step;
			n_channel = _n_channel;
			width_step = _width_step;

			data = (unsigned char *)_data;
		}
		PImage()
		{
			width =
				height =
				depth =
				n_channel =
				width_step = 0;

			data = 0;
		}

		PSize size() const
		{
			return PSize(width, height);
		}

		operator bool()
		{
			return 0 != data;
		}

	};

	inline int width(const PImage &img)
	{
		return img.width;
	}

	inline int height(const PImage &img)
	{
		return img.height;
	}

	inline PSize size(const PImage &img)
	{
		return img.size();
	}


	inline unsigned int n_channel(const PImage &img)
	{
		return img.n_channel;
	}

	inline unsigned int depth(const PImage &img)
	{
		return img.depth;
	}

	inline unsigned int width_step(const PImage &img)
	{
		return img.width_step;
	}

	inline unsigned char * scanline(const PImage &img,int y)
	{
		return img.data + y * (int)img.width_step;
	}

	inline unsigned char * scanline_bound(const PImage &img, int y)
	{
		if (y < 0) y = 0;
		if (y >= (int)img.height)y = img.height - 1;
		return scanline(img, y);
	}

	inline unsigned char * pick_pixel(const PImage &img, int x, int y)
	{
		return scanline(img, y) + x * (img.depth >> 3) * img.n_channel;
	}

	void copy(const PImage &dst, const PImage &src);

	void zero(const PImage &img);

	void zero(const PImage &img,unsigned int ch);

	void assert_range(const PImage &img, void * ptr);

	PImage clip_image(const PImage &img, int x_offset, int y_offset, int width, int height);

	template<typename DUnit>
	struct _COPY_CHANNEL_DFC
	{
		template<typename SUnit>
		DUnit operator () (SUnit v)
		{
			return (DUnit)v;
		}
	};

	template<typename DUnit,typename SUnit,typename C = _COPY_CHANNEL_DFC<DUnit> >
	void copy_channel(const PImage &dst, unsigned int dst_ch, const PImage &src, unsigned int src_ch, C && c = C())
	{
		assert(
			dst.width == src.width &&
			dst.height == src.height &&
			dst_ch < dst.n_channel &&
			src_ch < src.n_channel);

		for (int y = 0; y < dst.height; ++y)
		{
			DUnit * drow = (DUnit *)scanline(dst, y);
			SUnit * srow = (SUnit *)scanline(src, y);

			for (int x = 0; x < dst.width; ++x)
			{
				drow[dst_ch] = std::forward<C>(c)(srow[src_ch]);

				drow += dst.n_channel;
				srow += src.n_channel;
			}
		}
	}



#ifdef CV_MAJOR_VERSION

	inline PImage cv_use_roi(const IplImage & img)
	{
		if (!img->roi)return to_p_image(img);

		PImage ret;

		ret.width = img.roi->width;
		ret.height = img.roi->height;
		ret.n_channel = img.nChannels;
		ret.depth = img.depth;
		ret.width_step = img.widthStep;
		ret.data = (unsigned char *)(
			img->imageData +
			img->roi->yOffset * img->widthStep +
			img->roi->xOffset * (src->depth >> 3) * src->nChannels);

		return ret;
	}

	inline PImage cv_use_roi(const IplImage * img)
	{
		return cv_use_roi(*img);
	}

	inline cv::Mat to_cv_mat(const PImage &img)
	{
		int type;
		if (8 == img.depth)
		{
			type = CV_MAKETYPE(CV_8U, img.n_channel);
		}
		if (16 == depth(img))
		{
			type = CV_MAKETYPE(CV_16U, img.n_channel);
		}
		if (32 == depth(img))
		{
			type = CV_MAKETYPE(CV_32F, img.n_channel);
		}
		cv::Mat ret(img.height, img.width, type, img.date, img.width_step);
		return ret;
	}

	inline CvMat to_cv_cvmat(const PImage &img)
	{
		return to_cv_mat(img);
	}

	inline IplImage to_cv_image(const PImage &img)
	{
		return to_cv_mat(img);;
	}

#endif



	struct Image : public PImage
	{

		Image(
			unsigned int _width,
			unsigned int _height,
			unsigned int _n_channel,
			unsigned int _depth,
			unsigned int align = 4)
		{
			buf = 0;
			construct(_width, _height, _n_channel, _depth, align);
		}

		Image(
			PSize size,
			unsigned int _n_channel,
			unsigned int _depth,
			unsigned int align = 4)
		{
			buf = 0;
			construct(size, _n_channel, _depth, align);
		}
		Image()
		{
			buf = 0;
		}
		~Image()
		{
			destruct();
		}

		Image(const Image&) = delete;
		Image &operator =(const Image&) = delete;
		Image(Image &&other);
		Image & operator =(Image&& other);

		void construct(
			unsigned int _width,
			unsigned int _height,
			unsigned int _n_channel,
			unsigned int _depth,
			unsigned int align = 4);

		inline void construct(
			PSize size,
			unsigned int _n_channel,
			unsigned int _depth,
			unsigned int align = 4)
		{
			construct(size.width, size.height,
				_n_channel, _depth, align);
		}

		void destruct();

	private:


		unsigned char * buf;
	};

}

#endif
