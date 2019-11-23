#include <vector>
#include <iostream>
#include <tiffio.h>
#include <jpeglib.h>
#include <string>
#include <cmath>
#include <bear/image.h>
#include <bear/dynamic_image.h>
#include <bear/functor.h>
#include <boost/filesystem.hpp>
#include "down_semple.h"

using namespace std;
using namespace bear;
using namespace boost::filesystem;


using image_t = image<unsigned char, 3>;

image_t make_preview(
	path main_path,
	const_string_ptr file,
	const_array_ptr<size_t> x_grid,
	const_array_ptr<size_t> y_grid,
	size_t redundance,
	size_t sub_divide,
	size_t dst_width,  
	size_t dst_height)
{
	if (!exists(main_path))
	{
		return image_t();
	}

	auto full_width = x_grid.back();
	auto full_height = y_grid.back();
	auto full_x = x_grid.size() - 1;
	auto full_y = y_grid.size() - 1;

	tensor<image_t, 2> images(full_y, full_x);
	auto found = false;

	for (int y = 0; y < full_y; y++)
	{
		size_t y_left = y_grid[y];
		size_t y_right = y_grid[y + 1];
		size_t oh = y_right - y_left;
		if (y != 0)
		{
			oh += redundance;
		}
		if (y != full_y - 1)
		{
			oh += redundance;
		}

		for (int x = 0; x < full_x; x++)
		{
			auto base_path = main_path / ("_" + to_string(x) + "_" + to_string(y));
			auto file_path = base_path / string(file);
			size_t x_left = x_grid[x];
			size_t x_right = x_grid[x + 1];
			size_t ow = x_right - x_left;
			if (x != 0)
			{
				ow += redundance;
			}
			if (x != full_x - 1)
			{
				ow += redundance;
			}

			try
			{
				if (!exists(file_path))
				{
					throw bear_exception(exception_type::other_error, "");
				}
				TIFF *tiff = TIFFOpen(file_path.string().c_str(), "r");
				if (tiff == NULL)
				{
					throw bear_exception(exception_type::other_error, "unknown file io failed!");
				}
				defer df([=]() {
					TIFFClose(tiff);
				});

				uint32 w, h;
				uint16 cn, depth;

				TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
				TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
				TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &cn);
				TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &depth);

				if (cn < 3 || (depth != 8 && depth != 16))
				{
					throw bear_exception(
						exception_type::other_error,
						"unsupported format:",
						to_string(cn),
						" ",
						to_string(depth));
				}

				if (w != ow || h != oh)
				{
					throw bear_exception(
						exception_type::other_error,
						"wrong image size:",
						to_string(x),
						" ",
						to_string(y),
						" ",
						to_string(w),
						" ",
						to_string(h));
				}

				image_t oimg(w, h);
				dynamic_image dimg;
				dynamic_image_ptr pimg;

				if (cn != 3 || depth != 8)
				{
					dimg = dynamic_image(w, h, cn, image_unsigned_int_type, depth / 8);
					pimg = dimg;
				}
				else
				{
					pimg = oimg;
				}


				for (int m = 0; m < (int)h; ++m)
				{
					int suc = TIFFReadScanline(tiff, scanline(pimg, m).data(), m);
					if (suc != 1)
					{
						throw bear_exception(
							exception_type::other_error,
							"read image failed:",
							to_string(x),
							" ",
							to_string(y),
							" ",
							to_string(m));
					}
				}

				if (cn != 3 || depth != 8)
				{
					switch (depth)
					{
					case 8:
					{
						auto p_img = tensor_ptr<unsigned char, 3>(pimg);

						zip_to<2>([](auto d, auto s)
						{
							d[0] = s[0];
							d[1] = s[1];
							d[2] = s[2];
						}, oimg, p_img);
					}
					break;
					case 16:
					{
						auto p_img = tensor_ptr<unsigned short, 3>(pimg);

						zip_to<2>([](auto d, auto s)
						{
							d[0] = s[0] >> 8;
							d[1] = s[1] >> 8;
							d[2] = s[2] >> 8;
						}, oimg, p_img);
					}
					break;
					}
				}

				images[y][x] = move(oimg);
				found = true;


			}
			catch (bear_exception e)
			{
				if (!e.what().empty())
				{
					cout << "error:" << e.what() << endl;
				}
				vector<size_t> sub_x_grid;
				for (int i = 0; i < sub_divide + 1; ++i)
				{
					sub_x_grid.push_back(i * ow / sub_divide);
				}

				vector<size_t> sub_y_grid;
				for (int i = 0; i < sub_divide + 1; ++i)
				{
					sub_y_grid.push_back(i * oh / sub_divide);
				}

				images[y][x] = make_preview(
					base_path,
					file,
					sub_x_grid,
					sub_y_grid,
					redundance,
					sub_divide,
					ow,
					oh);

				if (!images[y][x].empty())
				{
					found = true;
				}
			}
		}
	}


	if (!found)
	{
		return image_t();
	}

	image_t dst(dst_width, dst_height);

	for (size_t y = 0; y < height(images); ++y)
	{
		auto oyl = y_grid[y];
		auto oyr = y_grid[y + 1];

		auto yofs = 0;
		auto ylmt = 0;
		if (y != 0)
		{
			yofs = 2;
		}
		if (y != height(images) - 1)
		{
			ylmt = 2;
		}

		auto ryl = oyl - yofs;
		auto ryr = oyr + ylmt;

		auto dyl = oyl * dst_height / y_grid.back();
		auto dyr = oyr * dst_height / y_grid.back();


		for (size_t x = 0; x < width(images); ++x)
		{
			if (images[y][x].empty())
			{
				continue;
			}
			auto oxl = x_grid[x];
			auto oxr = x_grid[x + 1];

			auto xofs = 0;
			auto xlmt = 0;
			if (x != 0)
			{
				xofs = 2;
			}
			if (x != width(images) - 1)
			{
				xlmt = 2;
			}

			auto rxl = oxl - xofs;
			auto rxr = oxr + xlmt;

			auto dxl = oxl * dst_width / x_grid.back();
			auto dxr = oxr * dst_width / x_grid.back();

			auto di = clip_image(dst, image_rectangle(dxl, dyl, dxr - dxl, dyr - dyl));

			if (dst_width == x_grid.back() && dst_height == y_grid.back())
			{
				auto si = clip_image(images[y][x], image_rectangle(xofs, yofs, dxr - dxl, dyr - dyl));
				copy(di, si);
				continue;
			}

			float x_fac = (float)dst_width / x_grid.back();
			float y_fac = (float)dst_height / y_grid.back();

			float x_ofs = (dxl + 0.5f) / x_fac - rxl;
			float y_ofs = (dyl + 0.5f) / y_fac - ryl;

			down_semple(di, images[y][x], x_ofs, x_fac, y_ofs, y_fac);
		}
	}

	return dst;
}

void save_jpeg(const string &path, image_ptr<unsigned char, 3> img)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	int row_stride = 0;
	FILE* fp = NULL;
	JSAMPROW row_pointer[1];

	cinfo.err = jpeg_std_error(&jerr);

	jpeg_create_compress(&cinfo);
	fp = fopen(path.c_str(), "wb");
	if (fp == NULL)
	{
		throw bear_exception(exception_type::other_error, "create result failed");
	}
	defer df([=]() {
		fclose(fp);
	});
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width = int(width(img));
	cinfo.image_height = int(height(img));
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;//设置输入格式

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 80, 1);  // todo 1 == true
	jpeg_start_compress(&cinfo, TRUE);
	row_stride = int(width(img)) * cinfo.input_components;

	int y = 0;
	while (cinfo.next_scanline < cinfo.image_height)
	{
		row_pointer[0] = &img[y][0][0];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
		++y;
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
}

int main(int argc, char *argv[])
{
	try
	{
		if (argc < 9)
		{
			throw bear_exception(exception_type::other_error, "wrong argument!");
		}

		const_string_ptr _path = argv[1];
		const_string_ptr _file = argv[2];
		const_string_ptr _result_path = argv[3];
		const_string_ptr _x_grid = argv[4];
		const_string_ptr _y_grid = argv[5];
		const_string_ptr _redundance = argv[6];
		const_string_ptr _sub_divide = argv[7];
		const_string_ptr _max_size = argv[8];

		const_string_ptr path = _path;
		const_string_ptr result_path = _result_path;
		const_string_ptr file = _file;
		vector<size_t> x_grid = map_function(
			[](auto s) {
				return string_cast<size_t>(s);
			},
			split(_x_grid, '_'));
		vector<size_t> y_grid = map_function(
			[](auto s) {
				return string_cast<size_t>(s);
			},
			split(_y_grid, '_'));
		auto redundance = stoi(_redundance);
		auto sub_divide = stoi(_sub_divide);
		auto max_size = stoi(_max_size);

		auto ow = x_grid.back();
		auto oh = y_grid.back();
		auto dw = ow;
		auto dh = oh;

		if (ow > oh)
		{
			if (ow > max_size)
			{
				dw = max_size;
				dh = oh * dw / ow;
			}
		}
		else
		{
			if (oh > max_size)
			{
				dh = max_size;
				dw = ow * dh / oh;
			}
		}

		auto dst = make_preview(
			string(path),
			file,
			x_grid,
			y_grid,
			redundance,
			sub_divide,
			dw,
			dh);

		if (dst.empty())
		{
			return 0;
		}

		save_jpeg(result_path, dst);
	}
	catch (const exception &e)
	{
		cerr << e.what();
		return -1;
	}
	catch (const bear_exception &e)
	{
		cerr << e.what();
		return -1;
	}
	return 0;
}
