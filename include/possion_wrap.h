#ifndef _POSSION_WRAP_H
#define _POSSION_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif
	typedef void* imagePtr;
	typedef void* matImagePtr;
	void go_poisson_stiching(imagePtr dstPtr, matImagePtr srcPtr,unsigned int rd,unsigned int format);

	imagePtr newImagePtr(size_t width, size_t height, size_t channel_size, unsigned int eletype, size_t ele_size, char* data, size_t width_step);
	void freeImagePtr(imagePtr p);

	matImagePtr newMatImagePtr(unsigned int y_max);
	void matPushImage(matImagePtr m,unsigned int y, size_t width, size_t height, size_t channel_size, unsigned int eletype, size_t ele_size, char* data, size_t width_step);
	void freeMat(matImagePtr m);


#ifdef __cplusplus
}
#endif

#endif