#include <vector>
#include <iostream>
#include <string>
#include <bear/ptr_algorism.h>
#include <bear/functor.h>
#include <opencv2/opencv.hpp>

#include "../include/possion_solver.h"
#include "../include/possion_stiching.h"


using namespace std;
using namespace bear;

int main(int argc, char *argv[])
{
	try
	{
		image_debug = [](const_dynamic_image_ptr _img)
		{
			image<unsigned char, 3> img(_img.width(),_img.height());
			if (_img.elm_size() == 2)
			{
				zip_to<2>([](auto& d, auto s)
					{
						if (s == 0)
						{
							d[0] = 0;
							d[1] = 0;
							d[2] = 0;
						}
						else
						{
							d[0] = 255;
							d[1] = 255;
							d[2] = 255;
						}
					}, img, const_image_ptr<unsigned short, 1>(_img));
			}
			else if (_img.channel_size() == 1)
			{
				zip_to<2>([](auto& d, auto s)
					{
						d[0] = s;
						d[1] = s;
						d[2] = s;
					}, img, const_image_ptr<unsigned char,1>(_img));
			}
			else
			{
				copy(img, const_image_ptr<unsigned char, 3>(_img));
			}
			cv::imshow("debug", cv::Mat(dynamic_image_ptr(img)));
			cv::waitKey();
		};

		auto _base = cv::imread("base.png");
		auto _patch = cv::imread("patch.png");
		auto _mask = cv::imread("mask.png");

		auto mask = image_cast<image_ptr<unsigned char, 3>>(_mask);

		mask.to_tensor_3d().for_each([](auto& p)
			{
				if (p > 128)
				{
					p = 0;
				}
				else
				{
					p = 255;
				}
			});

		PStichingParam param;

		param.constrain = PossionMaskConstrain;
		param.iteration_time = 100;

		poisson_stiching(_mask, _base, _patch, _mask, param);

		cv::imshow("111", _mask);
		cv::waitKey();
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