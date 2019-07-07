#ifndef _POSSION_WRAP_H
#define _POSSION_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif
	typedef void* imagePtr;
	typedef void* matImagePtr; 
	void go_poisson_stiching(imagePtr dstPtr, matImagePtr srcPtr,unsigned int rd,unsigned int format, unsigned int mode);

	imagePtr newImagePtr(unsigned long long width, unsigned long long height, unsigned long long channel_size, unsigned int eletype, unsigned long long ele_size, char* data, unsigned long long width_step);
	void freeImagePtr(imagePtr p);

	matImagePtr newMatImagePtr(unsigned int y_max);
	unsigned int matPushImage(matImagePtr m,unsigned int y, unsigned long long width, unsigned long long height, unsigned long long channel_size, unsigned int eletype, unsigned long long ele_size, char* data, unsigned long long width_step);
	void freeMat(matImagePtr m);


#ifdef __cplusplus
}
#endif

#endif