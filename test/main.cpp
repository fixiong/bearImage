#include <opencv2/opencv.hpp>
#include <types.h>
#include <image.h>
#include <possion_stiching.h>

using namespace cv;
using namespace std;
using namespace bear;

int main(){



    printf("Main Entry\n");

    int tileX=3,tileY=4;

    int w=320;
    int h=180;
    Mat dst(Size(tileX*w,tileY*h),CV_8UC3,Scalar(0));
    
    Image msk(w * 3, h * 4, 1, 8);
    for(int i=0;i<tileX;i++){
        for(int j=0;j<tileY;j++){
            char path[128];
            sprintf(path,"/Users/john/SandBox/done/12a06a1bb099fc8a770071471037a6a6.zip_t_%d_%d_res.png",i+1,j);
            Mat src=imread(path);

            imshow("src",src);
            waitKey(200);
            
            PImage dm = clip_image(dst, i * w, (3 - j) * h, w, h);
            PImage sm = src;
            copy(dm, sm);

            PImage ms = clip_image(msk, i * w, (3 - j) * h, w, h);

            if ((i + j * 3) & 1)
            {
                zero(ms);
            }
            else
            {
                fill(ms, 0, (unsigned char)255);
            }
            
        }
    }

    imshow("msk",to_cv_mat(msk));
    poisson_stiching_merged(dst, dst, msk, F_BGR, 5);

    imshow("dst",dst);
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