#include <vector>
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <string>
#include <bear/ptr_algorism.h>

#include "../include/possion_stiching.h"


using namespace cv;
using namespace std;
using namespace bear;

int main(){

	try
	{
		set_debug_callback([](dynamic_image_ptr _img, const_string_ptr) {
			Mat img = _img;
			imshow("debug", img);
			waitKey();
		});

		printf("Main Entry\n");

		int tileX = 4, tileY = 5;
		int rd = 2;


		vector<Mat> border(tileY);

		for (int y = 0; y < tileY; y++)
		{
			string name = "C:\\work\\images\\_r_";
			name += to_string(y);
			name += "\\result.tiff";
			border[y] = imread(name.c_str(), IMREAD_ANYDEPTH | IMREAD_COLOR);
		}

		auto ab = map_function([](const Mat & m)
		{
			return const_dynamic_image_ptr(m);
		}, border);

		vector<vector<Mat>> src(tileY, vector<Mat>(tileX));

		for (int x = 0; x < tileX; x++)
		{
			for (int y = 0; y < tileY; y++)
			{
				string name = "C:\\work\\images\\_";
				name += to_string(x);
				name += '_';
				name += to_string(y);
				name += "\\result.tiff";
				try {
					src[y][x] = imread(name.c_str(), IMREAD_ANYDEPTH | IMREAD_COLOR);
				}
				catch (...) {}
			}
		}
		int tw = rd * 2;
		int th = rd * 2;

		for (int x = 0; x < tileX; x++)
		{
			tw += src[0][x].cols - rd * 2;
		}

		for (int y = 0; y < tileY; y++)
		{
			th += src[y][0].rows - rd * 2;
		}

		Mat dst(Size(tw, th), CV_8UC3, Scalar(0));

		PStichingParam param;
		param.iteration_time = 100;


		Mat border_mat(Size(2, th), CV_8UC3, Scalar(0));
		param.panorama_border = border_mat;
		make_panorama_border(border_mat, ab);
		param.constrain = PossionPanoramaBorderConstrain;
		vector<size_t> x_grid(tileX + 1);
		vector<size_t> y_grid(tileY + 1);
		calculate_grid(src, x_grid, y_grid, rd);
		poisson_stiching(dst, src, x_grid, y_grid, rd, param);

		cout << endl;
		for_each(to_ptr(x_grid), [](auto s) 
		{
			cout << s << "_";
		});
		cout << endl;
		for_each(to_ptr(y_grid), [](auto s)
		{
			cout << s << "_";
		});

		//vector<image_point> error;

		//poisson_stiching_check(error, vector<image_point>(), src, rd, F_BGR, 8, 30);

		//for (int i = 0; i < (int)error.size(); ++i)
		//{
		//	cout << "error block: (" << error[i].x << "," << error[i].y << ")" << endl;
		//}

		imshow("aaa", dst);
		imwrite("c:\\work\\result.png", dst);
		waitKey();
	}
	catch (const cv::Exception &e)
	{
		cout << e.err;
	}
	catch (const bear_exception & e)
	{
		cout << e.what() << endl;
	}

	return 0;
}