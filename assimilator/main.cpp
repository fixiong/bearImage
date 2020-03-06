#include <vector>
#include <iostream>
#include <tiffio.h>
#include <string>
#include <bear/ptr_algorism.h>
#include <bear/functor.h>
#include <experimental/filesystem>

#include "../include/possion_stiching.h"

enum Mode
{
	ModeNormal,
	ModePanorama,
};

struct GlobalFormat {
	bool inited;
	uint16 photoMetric, planarConfig;
	uint16 cn, depth;
};

using namespace std;
using namespace bear;
using namespace std::experimental::filesystem::v1;

dynamic_image assimulate(
	path main_path,
	const_string_ptr file,
	const_array_ptr<size_t> x_grid,
	const_array_ptr<size_t> y_grid,
	size_t redundance,
	size_t sub_divide,
	GlobalFormat &format,
	Mode mode,
	size_t border_size,
	size_t border_width,
	bool complete)
{
	if (!exists(main_path))
	{
		return dynamic_image();
	}

	auto full_width = x_grid.back();
	auto full_height = y_grid.back();
	auto full_x = x_grid.size() - 1;
	auto full_y = y_grid.size() - 1;

	tensor<dynamic_image, 2> images(full_y, full_x);
	vector<dynamic_image> border(border_size);

	auto mx = full_x;
	if (mode == ModePanorama)
	{
		mx += 1;
	}

	auto found = false;
	auto lost = false;

	for (int y = 0; y < full_y; y++)
	{
		size_t bh = y_grid[y + 1] - y_grid[y];
		size_t rh = bh;
		if (y != 0)
		{
			rh += redundance;
		}
		if (y != full_y - 1)
		{
			rh += redundance;
		}

		for (int x = 0; x < mx; x++)
		{
			size_t bw = x_grid[x + 1] - x_grid[x];
			size_t rw = bw;
			path base_path;
			if (x == full_x)
			{
				base_path = main_path / ( "_r_" + to_string(y));
				bw = rw = border_width;
				rh = bh;
			}
			else
			{
				base_path = main_path / ("_" + to_string(x) + "_" + to_string(y));
				if (x != 0)
				{
					rw += redundance;
				}
				if (x != full_x - 1)
				{
					rw += redundance;
				}
			}

			auto file_path = base_path / string(file);

			dynamic_image img;

			
			if (!exists(file_path))
			{
				vector<size_t> sub_x_grid;
				for (int i = 0; i < sub_divide + 1; ++i)
				{
					sub_x_grid.push_back(i * rw / sub_divide);
				}

				vector<size_t> sub_y_grid;
				for (int i = 0; i < sub_divide + 1; ++i)
				{
					sub_y_grid.push_back(i * rh / sub_divide);
				}

				img = assimulate(
					base_path,
					file,
					sub_x_grid,
					sub_y_grid,
					redundance,
					sub_divide,
					format,
					ModeNormal,
					0,
					0,
					true);

			}
			else
			{
				TIFF *tiff = TIFFOpen(file_path.string().c_str(), "r");
				if (tiff == NULL)
				{
					cout << "read file failed!";
					lost = true;
					continue;
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

				if (w != rw || h != rh)
				{
					cout << "wrong image size:"
						<< x << ' '
						<< y << ' '
						<< w << ' '
						<< h << endl;
					lost = true;
					continue;
				}

				if (format.inited)
				{
					if (cn != format.cn || depth != format.depth)
					{
						cout << "wrong image format:"
							<< x << " "
							<< y << " "
							<< cn << " "
							<< depth << " "
							<< format.cn << " "
							<< format.depth << endl;
						lost = true;
						continue;
					}
				}
				else
				{
					format.inited = true;
					TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &format.photoMetric);
					TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &format.planarConfig);
					format.cn = cn;
					format.depth = depth;
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

				img = dynamic_image(w, h, cn, image_unsigned_int_type, depth / 8);

				auto read_error = false;
				for (int m = 0; m < (int)h; ++m)
				{
					int suc = TIFFReadScanline(tiff, scanline(img, m).data(), m);
					if (suc != 1)
					{
						read_error = true;
						cout << "read image failed:"
							<< x << " "
							<< y << " "
							<< m << endl;
						break;
					}
				}
				if (read_error){
					lost = true;
					continue;
				}
			}

			if (img.empty())
			{
				lost = true;
				continue;
			}

			found = true;
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

	if (!found || (complete && lost))
	{
		return dynamic_image();
	}

	dynamic_image dst(full_width, full_height, format.cn, image_unsigned_int_type, format.depth / 8);

	if (mode == ModePanorama)
	{
		PStichingParam param;
		param.iteration_time = 100;
		param.constrain = PossionPanoramaBorderConstrain;
		dynamic_image bd(2, full_height, format.cn, image_unsigned_int_type, format.depth / 8);
		param.panorama_border = bd;
		auto bdpt = map_function([](const dynamic_image &img) {
			return const_dynamic_image_ptr(img);
		}, border);

		make_panorama_border(bd, bdpt);
		poisson_stiching(dst, images, x_grid, y_grid, redundance, param);
	}
	else
	{
		PStichingParam param;
		param.iteration_time = 100;
		param.constrain = PossionNoConstrain;

		poisson_stiching(dst, images, x_grid, y_grid, redundance, param);
	}

	return dst;
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
		const_string_ptr _mode = argv[8];
		const_string_ptr _border_size;
		const_string_ptr _border_width;
		if (argc > 10)
		{
			_border_size = argv[9];
			_border_width = argv[10];
		}

		auto path = _path;
		auto result_path = _result_path;
		auto file = _file;
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
		auto redundance = stoi(_redundance);
		auto sub_divide = stoi(_sub_divide);
		Mode mode;
		size_t border_size = 0;
		size_t border_width = 0;
		if (_mode == "normal")
		{
			mode = ModeNormal;
		}
		else if (_mode == "panorama")
		{
			if (_border_size.empty() || _border_width.empty())
			{
				throw bear_exception(exception_type::other_error, "wrong panamorma argument!");
			}
			mode = ModePanorama;
			border_size = stoi(_border_size);
			border_width = stoi(_border_width);
		}
		else
		{
			throw bear_exception(exception_type::other_error, "wrong mode!");
		}

		GlobalFormat format = { 0 };

		auto img = assimulate(
			string(path) + "/",
			file,
			x_grid,
			y_grid,
			redundance,
			sub_divide,
			format,
			mode,
			border_size,
			border_width,
			false);

		if (img.empty())
		{
			throw bear_exception(exception_type::other_error, "file not found!");
		}

		string save_path = result_path;
		save_path += "/";
		save_path += file;
		//cout << save_path;
		TIFF *out = TIFFOpen(save_path.c_str(), "w");
		if (out == NULL)
		{
			throw bear_exception(exception_type::other_error, "creat result failed!");
		}
		defer df([=]() {
			TIFFClose(out);
		});

		uint16 compression = COMPRESSION_LZW; //
		TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width(img));
		TIFFSetField(out, TIFFTAG_IMAGELENGTH, height(img));
		TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
		TIFFSetField(out, TIFFTAG_COMPRESSION, compression);
		TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, format.depth);
		TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, format.cn);
		TIFFSetField(out, TIFFTAG_PHOTOMETRIC, format.photoMetric);
		TIFFSetField(out, TIFFTAG_PLANARCONFIG, format.planarConfig);

		for (size_t y = 0; y < height(img); y++)
		{
			int suc = TIFFWriteScanline(out, scanline(img, y).data(), (unsigned int)y);
			if (suc != 1)
			{
				throw bear_exception(exception_type::other_error, "save result failed!");
			}
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