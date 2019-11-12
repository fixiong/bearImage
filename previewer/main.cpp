#include <vector>
#include <iostream>
#include <tiffio.h>
#include <string>
#include <bear/image.h>
#include <bear/dynamic_image.h>
#include <bear/functor.h>

using namespace std;
using namespace bear;

int main(int argc, char *argv[])
{
	try
	{
		if (argc < 10)
		{
			throw bear_exception(exception_type::other_error, "wrong argument!");
		}

		const_string_ptr _path = argv[1];
		const_string_ptr _file = argv[2];
		const_string_ptr _result_path = argv[3];
		const_string_ptr _full_x = argv[4];
		const_string_ptr _full_y = argv[5];
		const_string_ptr _x_grid = argv[6];
		const_string_ptr _y_grid = argv[7];
		const_string_ptr _redundance = argv[8];
		const_string_ptr _max_size = argv[9];

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
		auto max_size = stoi(_max_size);

		auto dw = x_grid.back();
		auto dh = y_grid.back();

		tensor<dynamic_image, 2> images(full_y, full_x);
		uint16 photoMetric, planarConfig;
		uint16 cn, depth;
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
				auto file_path = string(path) + "/_" + to_string(x) + "_" + to_string(y) + "/" + string(file);
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

				TIFF *tiff = TIFFOpen(file_path.c_str(), "r");
				if (tiff == NULL)
				{
					continue;
				}
				defer df([=]() {
					TIFFClose(tiff);
				});

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

				images[y][x] = move(img);
			}
		}

		dynamic_image dst(dw, dh, 3, image_unsigned_int_type, 1);

		if (!found)
		{

		}

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