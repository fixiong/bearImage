#include <vector>
#include <iostream>
#include <tiffio.h>
#include <jpeglib.h>
#include <string>
#include <bear/image.h>
#include <bear/dynamic_image.h>
#include <bear/functor.h>

using namespace std;
using namespace bear;


using image_t = image<unsigned char, 3>;

image_t make_preview(
	string path,
	const_string_ptr file,
	size_t full_x,
	size_t full_y,
	const_array_ptr<size_t> x_grid,
	const_array_ptr<size_t> y_grid,
	size_t redundance,
	size_t sub_divide,
	size_t dst_width,
	size_t dst_height)
{
	auto dw = x_grid.back();
	auto dh = y_grid.back();

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
			auto base_path = path + "/_" + to_string(x) + "_" + to_string(y);
			auto file_path = base_path + "/" + string(file);
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
					ow,
					oh,
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
			else
			{
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
					continue;
				}

				if (w != ow || h != oh)
				{
					continue;
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
		}
	}


	if (!found)
	{
		return image_t();
	}

	image_t dst(dst_width, dst_height);

	size_t current_y_grid = 0;
	for (size_t y = 0; y < height(dst); ++y)
	{
		size_t sy = y * dh / height(dst);
		while (y_grid[current_y_grid + 1] <= sy)
		{
			++current_y_grid;
		}
		sy -= y_grid[current_y_grid];
		if (current_y_grid != 0)
		{
			sy += redundance;
		}

		auto drow = dst[y];
		auto xrows = map_function([=](const image_t &imgs)
		{
			return imgs[sy];
		}, images[current_y_grid]);

		size_t current_x_grid = 0;
		for (size_t x = 0; x < width(dst); ++x)
		{
			size_t sx = x * dw / width(dst);
			while (x_grid[current_x_grid + 1] <= sx)
			{
				++current_x_grid;
			}
			sx -= x_grid[current_x_grid];
			if (current_x_grid != 0)
			{
				sx += redundance;
			}

			auto dp = drow[x];
			auto ds = xrows[current_x_grid][sx];

			dp[0] = ds[0];
			dp[1] = ds[1];
			dp[2] = ds[2];
		}
	}


	return dst;
}

//===================================================================================
//function:       jpegѹ��
//input:          1:���ɵ��ļ���,2:bmp��ָ��,3:λͼ���,4:λͼ�߶�,5:��ɫ���
//return:         int
//description:    bmp�����ظ�ʽΪ(RGB)
//===================================================================================
int savejpeg(char *filename, unsigned char *bits, int width, int height, int depth)
{
	struct jpeg_compress_struct jcinfo;  //����jpegѹ������
	struct jpeg_error_mgr jerr;
	FILE * outfile;                 //target file 
	JSAMPROW row_pointer[1];        //pointer to JSAMPLE row[s] һ��λͼ 
	int     row_stride;             //ÿһ�е��ֽ��� 
	jcinfo.err = jpeg_std_error(&jerr);   //ָ����������
	jpeg_create_compress(&jcinfo);      //��ʼ��jpegѹ������

	//ָ��ѹ�����ͼ������ŵ�Ŀ���ļ���ע�⣬Ŀ���ļ�Ӧ�Զ�����ģʽ��
	if ((outfile = fopen(filename, "wb")) == NULL)
	{
		fprintf(stderr, "can't open %s/n", filename);
		return -1;
	}
	jpeg_stdio_dest(&jcinfo, outfile);   //ָ��ѹ�����ͼ������ŵ�Ŀ���ļ�
	jcinfo.image_width = width;      // Ϊͼ�Ŀ�͸ߣ���λΪ����
	jcinfo.image_height = height;
	jcinfo.input_components = 3;         // �ڴ�Ϊ3,��ʾ��ɫλͼ�� ����ǻҶ�ͼ����Ϊ1
	jcinfo.in_color_space = JCS_RGB;         //JCS_GRAYSCALE��ʾ�Ҷ�ͼ��JCS_RGB��ʾ��ɫͼ��
	/*
	��Ҫע����ǣ�jpeg_set_defaults����һ��Ҫ�����ú�ͼ����ߡ�ɫ��ͨ������ɫ�ʿռ��ĸ���������ܵ��ã���Ϊ�������Ҫ�õ����ĸ�ֵ������jpeg_set_defaults������jpeglib �����Ĭ�ϵ����ö�ͼ�����ѹ��
	�����Ҫ�ı����ã���ѹ��������������������󣬿��Ե����������ú�������jpeg_set_quality��������ʵͼ��ѹ��ʱ�кö������������
	���󲿷����Ƕ��ò������ã�ֻ�����jpeg_set_defaults����ֵΪĬ��ֵ����
	*/
	jpeg_set_defaults(&jcinfo);
	jpeg_set_quality(&jcinfo, 60, TRUE);//limit to baseline-JPEG values
	/*
	���ȵ���jpeg_start_compress��Ȼ����Զ�ÿһ�н���ѹ����Ҳ���Զ������н���ѹ�����������Զ�������ͼ�����һ��ѹ����ѹ����ɺ󣬼ǵ�Ҫ����jpeg_finish_compress����
	*/
	jpeg_start_compress(&jcinfo, TRUE);

	row_stride = width * depth; // JSAMPLEs per row in image_buffer(���������ͼ����Ҫ����3) 
	//��ÿһ�н���ѹ��
	while (jcinfo.next_scanline < jcinfo.image_height)
	{
		//�����������޸ģ�����jpg�ļ���ͼ���ǵ��ģ����Ը���һ�¶���˳��
		//����ԭ���룺
		//row_pointer[0] = & bits[jcinfo.next_scanline * row_stride];
		row_pointer[0] = &bits[(jcinfo.image_height - jcinfo.next_scanline - 1) * row_stride];
		(void)jpeg_write_scanlines(&jcinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&jcinfo);
	fclose(outfile);
	jpeg_destroy_compress(&jcinfo);
	return 0;
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
	cinfo.in_color_space = JCS_RGB;//���������ʽ

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 50, 1);  // todo 1 == true
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
		if (argc < 11)
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
		const_string_ptr _sub_divide = argv[9];
		const_string_ptr _max_size = argv[10];

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
			path,
			file,
			full_x,
			full_y,
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
