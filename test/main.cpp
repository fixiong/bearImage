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

	set_debug_callback([](dynamic_image_ptr _img, const_string_ptr) {
		Mat img = _img;
		imshow("debug", img);
		waitKey();
	});

	printf("Main Entry\n");

	int tileX = 5, tileY = 5;
	int rd = 2;


	vector<Mat> border(tileY);

	for (int y = 0; y < tileY; y++)
	{
		string name = "C:\\work\\pictures\\4204\\_r_";
		name += to_string(y);
		name += "\\p\\result.tiff";
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
			string name = "C:\\work\\pictures\\4204\\_";
			name += to_string(x);
			name += '_';
			name += to_string(y);
			name += "\\p\\result.tiff";
			src[y][x] = imread(name.c_str(), IMREAD_ANYDEPTH | IMREAD_COLOR);
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


	try
	{
		Mat border(Size(2, th), CV_8UC3, Scalar(0));
		param.panorama_border = border;
		make_panorama_border(border, ab);
		param.constrain = PossionPanoramaBorderConstrain;
		poisson_stiching(dst, src, rd, param);
	}
	catch (const bear_exception & e)
	{
		cout << e.what() << endl;
		exit(0);
	}

	//vector<image_point> error;

	//poisson_stiching_check(error, vector<image_point>(), src, rd, F_BGR, 8, 30);

	//for (int i = 0; i < (int)error.size(); ++i)
	//{
	//	cout << "error block: (" << error[i].x << "," << error[i].y << ")" << endl;
	//}

	imshow("aaa", dst);
	imwrite("c:\\work\\result.png", dst);
	waitKey();



    return 0;
}