#include <iostream>
#include <tiffio.h>
#include <string>
#include <bear/ptr_algorism.h>
#include <bear/functor.h>

using namespace std;
using namespace bear;

int main(int argc, char *argv[])
{
	try
	{
		if (argc < 4)
		{
			throw bear_exception(exception_type::other_error, "wrong argument!");
		}

		const_string_ptr _path = argv[1];
		const_string_ptr _width = argv[2];
		const_string_ptr _height = argv[3];

		const_string_ptr path = _path;
		auto width = stoi(_width);
		auto height = stoi(_height);
		vector<const_string_ptr> files(argv + 4, argv + argc);

		for_each(files, [=](auto filename) {
			string full_path = path;
			full_path += '/';
			full_path += filename;

			TIFF* tiff = TIFFOpen(full_path.c_str(), "r");
			if (tiff == NULL)
			{
				throw bear_exception(exception_type::other_error, "open ", full_path," failed");
			}
			defer df([=]() {
				TIFFClose(tiff);
				});

			uint32 w, h;

			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);

			if (w != width || h != height) {
				throw bear_exception(
					exception_type::other_error,
					"wrong size:",
					full_path,
					" expect:",
					to_string(width),
					" ",
					to_string(height),
					" has:",
					to_string(w),
					" ",
					to_string(h));
			}

		});
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
