#ifndef _POSSION_WRAP_H
#define _POSSION_WRAP_H
#include "../include/bear/dynamic_image.h"

#ifdef __cplusplus
extern "C" {
#endif
	typedef void* imagePtr;
	typedef void* matImagePtr; 
	unsigned int assimilate(const char* workDirPath,const char* fileName, unsigned int tNumX, unsigned int tNumY, unsigned int mode, unsigned int param);
	void go_poisson_stiching(bear::dynamic_image_ptr dst, std::vector<std::vector<bear::dynamic_image_ptr>>  src, std::vector<std::vector<bear::dynamic_image_ptr>> border, unsigned int rd, unsigned int format, unsigned int mode);

	imagePtr newImagePtr(unsigned long long width, unsigned long long height, unsigned long long channel_size, unsigned int eletype, unsigned long long ele_size, char* data, unsigned long long width_step);
	void freeImagePtr(imagePtr p);
	
	matImagePtr newMatImagePtr(unsigned int y_max);
	unsigned int matPushImage(matImagePtr m,unsigned int y, unsigned long long width, unsigned long long height, unsigned long long channel_size, unsigned int eletype, unsigned long long ele_size, char* data, unsigned long long width_step);
	void freeMat(matImagePtr m);


#ifdef __cplusplus
}
#endif

#endif