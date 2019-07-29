#include "../include/possion_stiching.h"
#include "../include/possion_wrap.h"
#include "../include/bear/dynamic_image.h"
#include <iostream>
#include <string>

void go_poisson_stiching(imagePtr dstPtr, matImagePtr srcPtr, matImagePtr borderPtr, unsigned int rd, unsigned int format, unsigned int mode)
{
	try {
		std::vector<std::vector<bear::dynamic_image_ptr>>* src = (std::vector<std::vector<bear::dynamic_image_ptr>>*) srcPtr;
		std::vector<std::vector<bear::dynamic_image_ptr>>* border = (std::vector<std::vector<bear::dynamic_image_ptr>>*) borderPtr;
		bear::dynamic_image_ptr* dst = (bear::dynamic_image_ptr*) dstPtr;


		PStichingParam param;

		param.iteration_time = 100;
		if (mode)
		{
			param.constrain = PossionPanoramaBorderConstrain;

			if (!border || border->size() != 1 || border->at(0).empty()) {
				throw bear::bear_exception(bear::exception_type::other_error, "wrong panorama border!");
			}

			if (border->at(0)[0].elm_size() == 1)
			{
				bear::tensor<unsigned char, 3> bd(dst->height(), dst->width(), dst->channel_size());
				param.panorama_border = bd;

				PStichingVectorSrc bs(*border);
				make_panorama_border(bd, bs.src[0]);
				poisson_stiching(*dst, *src, rd, param);
			}
			else
			{
				bear::tensor<unsigned short, 3> bd(dst->height(), dst->width(), dst->channel_size());
				param.panorama_border = bd;

				PStichingVectorSrc bs(*border);
				make_panorama_border(bd, bs.src[0]);
				poisson_stiching(*dst, *src, rd, param);
			}
		}
		else
		{
			param.constrain = PossionNoConstrain;

			poisson_stiching(*dst, *src, rd, param);
		}
	}
	catch (bear::bear_exception e) {
		std::cout<<e.what()<<std::endl;
	}
}
imagePtr newImagePtr(unsigned long long width, unsigned long long height, unsigned long long channel_size, unsigned int eletype, unsigned long long ele_size, char* data, unsigned long long width_step)
{

	bear::dynamic_image_ptr* re;

	re = new bear::dynamic_image_ptr(width, height, channel_size, bear::data_type(eletype), ele_size, data, width_step);


	return (void*)re;
}
void freeImagePtr(imagePtr p)
{
	bear::dynamic_image_ptr* dip = (bear::dynamic_image_ptr*)p;
	delete dip;
}
matImagePtr newMatImagePtr(unsigned int y_max) 
{

	std::vector<std::vector<bear::dynamic_image_ptr>>* re = new std::vector<std::vector<bear::dynamic_image_ptr>>(y_max);

	//(*re)[0].push_back(image);

	return (void*)re;
}

unsigned int matPushImage(matImagePtr m,unsigned int y, unsigned long long width, unsigned long long height, unsigned long long channel_size, unsigned int eletype, unsigned long long ele_size, char* data, unsigned long long width_step)
{
	
	bear::dynamic_image_ptr image;
	try {
		image = bear::dynamic_image_ptr(width, height, channel_size, bear::data_type(eletype), ele_size, data, width_step);
	}
	catch (bear::bear_exception e) {
		std::cout << e.what() << std::endl;
		return 0;
	}
	std::vector<std::vector<bear::dynamic_image_ptr>>* v = (std::vector<std::vector<bear::dynamic_image_ptr>>*)m;
	(*v)[y].push_back(image);
	return 1;
}

void freeMat(matImagePtr m)
{
	std::vector<std::vector<bear::dynamic_image_ptr>>* mp = (std::vector<std::vector<bear::dynamic_image_ptr>>*)m;
	delete mp;
}
