#include <vector>
#include <iostream>
#include <tiffio.h>
#include <string>
#include <bear/ptr_algorism.h>

#include "../include/possion_stiching.h"

enum Mode
{
	ModeNormal,
	ModePananorma,
};

using namespace std;
using namespace bear;

int main(int _argc, char *_argv[])
{

	try
	{
		vector<const_string_ptr> arg;

		for (int i = 1; i < _argc; ++i)
		{
			auto a = const_string_ptr(_argv[i]);

			auto splt = a.split(' ');
			arg.push_back(splt.first);
			while (!splt.second.empty())
			{
				splt = splt.second.split(' ');
				arg.push_back(splt.first);
			}
		}

		if (arg.size() < 9)
		{
			throw bear_exception(exception_type::other_error, "wrong argument!");
		}

		auto _path = arg[0];
		auto _file = arg[1];
		auto _result_path = arg[2];
		auto _full_x = arg[3];
		auto _full_y = arg[4];
		auto _x_grid = arg[5];
		auto _y_grid = arg[6];
		auto _redundance = arg[7];
		auto _mode = arg[8];
		const_string_ptr _border_size;
		const_string_ptr _border_width;
		if (arg.size() > 10)
		{
			_border_size = arg[9];
			_border_width = arg[10];
		}

		auto path = _path;
		auto result_path = _result_path;
		auto file = _file;
		auto full_x = stoi(_full_x);
		auto full_y = stoi(_full_y);
		auto x_grid = map_function(
			[](auto s) {
				return string_cast<size_t>(s);
			},
			split(_x_grid, '_'));
		auto y_grid = map_function(
			[](auto s) {
				return string_cast<size_t>(s);
			},
			split(_y_grid, '_'));

		if (x_grid.size() != full_x + 1 || y_grid.size() != full_y + 1)
		{
			throw bear_exception(exception_type::other_error, "wrong grid size!");
		}
		auto redundance = stoi(_redundance);
		Mode mode;
		size_t border_size = 0;
		size_t border_width = 0;
		if (_mode == "normal")
		{
			mode = ModeNormal;
		}
		else if (_mode == "pananorma")
		{
			if (_border_size.empty() || _border_width.empty())
			{
				throw bear_exception(exception_type::other_error, "wrong panamorma argument!");
			}
			mode = ModePananorma;
			border_size = stoi(_border_size);
			border_width = stoi(_border_width);
		}
		else
		{
			throw bear_exception(exception_type::other_error, "wrong mode!");
		}

		auto dw = x_grid.back();
		auto dh = y_grid.back();

		tensor<dynamic_image, 2> images(full_y, full_x);
		uint16 photoMetric, planarConfig;
		uint16 cn, depth;
		auto found = false;

		vector<dynamic_image> border(border_size);

		auto mx = full_x;
		if (mode == ModePananorma)
		{
			mx += 1;
		}

		for (int y = 0; y < full_y; y++)
		{
			size_t y_left = y_grid[y];
			size_t y_right = y_grid[y + 1];
			size_t _oh = y_right - y_left;
			if (y != 0)
			{
				_oh += redundance;
			}
			if (y != full_y - 1)
			{
				_oh += redundance;
			}

			for (int x = 0; x < mx; x++)
			{
				size_t x_left;
				size_t x_right;
				size_t ow;
				size_t oh;
				string file_path;
				if (x == full_x)
				{
					file_path = string(path) + "/_r_" + to_string(y) + "/" + string(file);
					x_left = 0;
					x_right = border_width;
					ow = border_size;
					oh = y_right - y_left;
				}
				else
				{
					file_path = string(path) + "/_" + to_string(x) + "_" + to_string(y) + "/" + string(file);
					x_left = x_grid[x];
					x_right = x_grid[x + 1];
					ow = x_right - x_left;
					if (x != 0)
					{
						ow += redundance;
					}
					if (x != full_x - 1)
					{
						ow += redundance;
					}
					oh = _oh;
				}

				TIFF *tiff = TIFFOpen(file_path.c_str(), "r");
				if (tiff == NULL)
				{
					continue;
				}
				uint32 w, h;
				uint16 c_cn, c_depth;

				TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
				TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
				TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &c_cn);
				TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &c_depth);

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

				if (!found)
				{
					found = true;
					TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photoMetric);
					TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planarConfig);
					cn = c_cn;
					depth = c_depth;
					if (depth != 8 && depth != 16)
					{
						throw bear_exception(
							exception_type::other_error,
							"unsupported format:",
							to_string(cn),
							" ",
							to_string(depth));
					}
				}
				else if (cn != c_cn || depth != c_depth)
				{
					throw bear_exception(
						exception_type::other_error,
						"wrong image format:",
						to_string(x),
						" ",
						to_string(y),
						" ",
						to_string(cn),
						" ",
						to_string(depth),
						" ",
						to_string(c_cn),
						" ",
						to_string(c_depth));
				}

				dynamic_image img(w, h, cn, image_unsigned_int_type, depth / 8);

				for (int m = 0; m < (int)h; ++m)
				{
					int suc = TIFFReadScanline(tiff, scanline(img, m).data(), m);
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
				TIFFClose(tiff);

				if (x == full_x)
				{
					border[y] = move(img);
				}
				else
				{
					images[y][x] = move(img);
				}
			}
		}

		if (!found)
		{
			throw bear_exception(exception_type::other_error, "file not found!");
		}

		dynamic_image dst(dw, dh, cn, image_unsigned_int_type, depth / 8);

		if (mode == ModePananorma)
		{
			PStichingParam param;
			param.iteration_time = 100;
			param.constrain = PossionPanoramaBorderConstrain;

			if (depth == 8)
			{
				bear::tensor<unsigned char, 3> bd(dh, dw, cn);
				param.panorama_border = bd;

				make_panorama_border(bd, map_function(
											 [](const dynamic_image &img) {
												 return const_dynamic_image_ptr(img);
											 },
											 border));
				poisson_stiching(dst, images, x_grid, y_grid, redundance, param);
			}
			else
			{
				bear::tensor<unsigned short, 3> bd(dh, dw, cn);
				param.panorama_border = bd;

				make_panorama_border(bd, map_function(
											 [](const dynamic_image &img) {
												 return const_dynamic_image_ptr(img);
											 },
											 border));
				poisson_stiching(dst, images, x_grid, y_grid, redundance, param);
			}
		}
		else
		{
			PStichingParam param;
			param.constrain = PossionNoConstrain;

			poisson_stiching(dst, images, x_grid, y_grid, redundance, param);
		}

		string save_path = result_path;
		save_path += "/";
		save_path += file;
		TIFF *out = TIFFOpen(save_path.c_str(), "w");

		uint16 compression = COMPRESSION_LZW; //
		TIFFSetField(out, TIFFTAG_IMAGEWIDTH, dw);
		TIFFSetField(out, TIFFTAG_IMAGELENGTH, dh);
		TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
		TIFFSetField(out, TIFFTAG_COMPRESSION, compression);
		TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, depth);
		TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, cn);
		TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photoMetric);
		TIFFSetField(out, TIFFTAG_PLANARCONFIG, planarConfig);

		for (size_t m = 0; m < dh; m++)
		{
			int suc = TIFFWriteScanline(out, scanline(dst, m).data(), (unsigned int)m);
			if (suc != 1)
			{
				throw bear_exception(exception_type::other_error, "save result failed!");
			}
		}
		TIFFClose(out);
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