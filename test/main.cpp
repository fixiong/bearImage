//#include <opencv2/opencv.hpp>
#include <types.h>
#include <possion_stiching.h>
#include <vector>
#include <iostream>

//using namespace cv;
using namespace std;
using namespace bear;

int main(){



	//printf("Main Entry\n");

	//int tileX = 3, tileY = 4;
	//int rd = 3;

	//vector<vector<Mat>> src(tileY, vector<Mat>(tileX));


	//for (int x = 0; x < tileX; x++)
	//{
	//	for (int y = 0; y < tileY; y++)
	//	{
	//		char path[128];
	//		sprintf(path, "/Users/john/SandBox/done/12a06a1bb099fc8a770071471037a6a6.zip_t_%d_%d_res.png", x, y);
	//		src[y][x] = imread(path);

	//	}
	//}
	//int tw = rd * 2;
	//int th = rd * 2;

	//for (int x = 0; x < tileX; x++)
	//{
	//	tw += src[0][x].width - rd * 2;
	//}

	//for (int y = 0; y < tileY; y++)
	//{
	//	tw += src[y][0].height - rd * 2;
	//}

	//Mat dst(Size(tw, th), CV_8UC3, Scalar(0));

	//PStichingParam param;
	//param.iteration_time = 100;

	//poisson_stiching(dst, src, rd, F_BGR, 5, param);

	//vector<PPoint> error;

	//poisson_stiching_check(error, vector<PPoint>(), src, rd, F_BGR, 8, 30);

	//for (int i = 0; i < (int)error.size(); ++i)
	//{
	//	cout << "error block: (" << error[i].x << "," << error[i].y << ")" << endl;
	//}

	//imshow("dst", dst);
	//waitKey();



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