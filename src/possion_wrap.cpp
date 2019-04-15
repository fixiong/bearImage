#include "../include/possion_stiching.h"
#include "../include/possion_wrap.h"
#include "../include/bear/dynamic_image.h"

void go_poisson_stiching(imagePtr dstPtr, matImagePtr srcPtr, unsigned int rd, unsigned int format)
{
	auto src = (std::vector<std::vector<bear::dynamic_image_ptr>>*) srcPtr;
	auto dst = (bear::dynamic_image_ptr*) dstPtr;
	

	PStichingParam param;

	param.iteration_time = 100;
	param.constrain = PossionEstimateConstrain;

	poisson_stiching(*dst, *src, rd, format, param);

}
imagePtr newImagePtr(size_t width, size_t height, size_t channel_size, unsigned int eletype, size_t ele_size, char* data, size_t width_step)
{
	auto re = new bear::dynamic_image_ptr(width,height, channel_size, bear::data_type(eletype), ele_size, data, width_step);
	return (void*)re;
}
void freeImagePtr(imagePtr p)
{
	auto dip = (bear::dynamic_image_ptr*)p;
	delete dip;
}
matImagePtr newMatImagePtr(unsigned int y_max) 
{
	auto re = new std::vector<std::vector<bear::dynamic_image_ptr>>(y_max);
	return (void*)re;
}

void matPushImage(matImagePtr m,unsigned int y, size_t width, size_t height, size_t channel_size, unsigned int eletype, size_t ele_size, char* data, size_t width_step)
{
	bear::dynamic_image_ptr image(width, height, channel_size, bear::data_type(eletype), ele_size, data, width_step);
	auto v = (std::vector<std::vector<bear::dynamic_image_ptr>>*)m;
	(*v)[y].push_back(image);
}

void freeMat(matImagePtr m)
{
	auto mp = (std::vector<std::vector<bear::dynamic_image_ptr>>*)m;
	delete mp;
}