//#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <string>

#include "../include/possion_stiching.h"

using namespace cv;
using namespace std;
using namespace bear;

void show_debug(dynamic_image_ptr _img)
{
	Mat img = _img;
	imshow("dst", img);
	waitKey();
}

int main(){

	set_debug_callback(show_debug);

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
	param.iteration_time = 4;

	poisson_stiching(dst, src, rd, F_BGR, param);

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




// QImage img[4][3];

//  for (int y = 0; y < 4; ++y)
//  {
//   for (int x = 0; x < 3; ++x)
//   {
//    QString fn = QString("C:\\work\\pictures\\Archive\\") + (char)('1' + x) + '_' + (char)('0' + y) + "_res.png";
//    if (!img[y][x].load(fn)) //加载图像
//    {
//     QMessageBox::information(this,
//      tr("打开图像失败"),
//      tr("打开图像失败!"));
//     return;
//    }

//    if (img[y][x].format() != QImage::Format_RGBA8888)
//    {
//     img[y][x] = img[y][x].convertToFormat(QImage::Format_RGBA8888, Qt::ColorOnly);
//    }
//   }
//  }

//  int w = bear::width(img[0][0]);
//  int h = bear::height(img[0][0]);

//  QImage dst(QSize(w * 3,h * 4), QImage::Format_RGBA8888);

//  Image msk(w * 3, h * 4, 1, 8);

//  for (int y = 0; y < 4; ++y)
//  {
//   for (int x = 0; x < 3; ++x)
//   {
//    PImage dm = clip_image(dst, x * w, (3 - y) * h, w, h);
//    PImage sm = img[y][x];

//    copy(dm, sm);

//    PImage ms = clip_image(msk, x * w, (3 - y) * h, w, h);

//    if ((x + y * 3) & 1)
//    {
//     zero(ms);
//    }
//    else
//    {
//     fill(ms, 0, (unsigned char)255);
//    }
//   }
//  }

//  poisson_stiching_merged(dst, dst, msk, F_RGBA, 5);

//  ui.label->setPixmap(QPixmap::fromImage(dst));