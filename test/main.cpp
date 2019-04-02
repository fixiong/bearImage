#include <vector>
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <string>

#include "../include/possion_stiching.h"

using namespace cv;
using namespace std;
using namespace bear;

int main(){

	set_debug_callback([](dynamic_image_ptr _img, const_string_ptr) {
		Mat img = _img;
		imshow("dst", img);
		waitKey();
	});

	printf("Main Entry\n");

	int tileX = 2, tileY = 2;
	int rd = 2;

	vector<vector<Mat>> src(tileY, vector<Mat>(tileX));


	for (int x = 0; x < tileX; x++)
	{
		for (int y = 0; y < tileY; y++)
		{
			string name = "c:/work/images/img_";
			name += to_string(x);
			name += '_';
			name += to_string(y);
			name += ".png";
			src[y][x] = imread(name.c_str());

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
	param.constrain = PossionEstimateConstrain;

	//try
	//{
		poisson_stiching(dst, src, rd, F_BGR, param);
	//}
	//catch (const bear_exception & e)
	//{
	//	cout << e.what() << endl;
	//	exit(0);
	//}

	//vector<image_point> error;

	//poisson_stiching_check(error, vector<image_point>(), src, rd, F_BGR, 8, 30);

	//for (int i = 0; i < (int)error.size(); ++i)
	//{
	//	cout << "error block: (" << error[i].x << "," << error[i].y << ")" << endl;
	//}

	imshow("dst", dst);
	waitKey();



    return 0;
}