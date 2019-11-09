#include "../include/possion_stiching.h"
#include "../include/possion_wrap.h"
//#include "../include/bear/dynamic_image.h"
#include "tiffio.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;

void go_poisson_stiching(bear::dynamic_image_ptr* dst, std::vector<std::vector<bear::dynamic_image_ptr>>* src, std::vector<std::vector<bear::dynamic_image_ptr>>* border, unsigned int rd, unsigned int format, unsigned int mode)
{
	try {



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
				bear::tensor<unsigned char, 3> bd(dst->height(), 2, dst->channel_size());
				param.panorama_border = bd;

				PStichingVectorSrc bs(*border);
				make_panorama_border(bd, bs.src[0]);
				poisson_stiching(*dst, *src, rd, param);
			}
			else
			{
				bear::tensor<unsigned short, 3> bd(dst->height(), 2, dst->channel_size());
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

unsigned int assimilate(const char* workDirPath ,const char* fileName ,unsigned int tNumX, unsigned int tNumY, unsigned int mode, unsigned int param)
{
	uint16 cn, depth, photoMetric, planarConfig;
	uint32 width=0, height=0;
	vector<vector<bear::tensor<char, 3>>> mat(tNumY,vector<bear::tensor< char, 3>>(tNumX));
	vector<vector<bear::dynamic_image_ptr>> srcMat(tNumY);

	for (int idxResJ = 0; idxResJ < tNumY; idxResJ++)
	{
		for (int idxResI = 0; idxResI < tNumX; idxResI++)
		{
			string picFilePath = string(workDirPath) + "/_";
			picFilePath += to_string(idxResI) + "_" + to_string(idxResJ);
			picFilePath += "/channels/";
			picFilePath += fileName;

			TIFF* tiff = TIFFOpen(picFilePath.c_str(), "r");
			if (tiff == NULL){
				return 0;
			}
			uint32 w, h;
		
			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
			TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &cn);
			TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &depth);
			TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photoMetric);
			TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planarConfig);

			std::cout << "width\n" << w;
			std::cout << "height\n" << h;
			std::cout << "cn\n" << cn;
			std::cout << "depth\n" << depth;
			std::cout << "photoMetric\n" << photoMetric;
			std::cout << "planarConfig\n" << planarConfig;

			if (idxResI == 0) {
				height += h - 4;
			}
			if (idxResJ == 0) {
				width += w - 4;
			}
			if (depth > 8) {
				mat[idxResJ][idxResI] = bear::tensor<char, 3>(h, w, cn*2);
			}
			else {
				mat[idxResJ][idxResI] = bear::tensor<char, 3>(h, w, cn);
			}
	
			for (int m = 0; m < h; m++)
			{
				int suc = TIFFReadScanline(tiff, &mat[idxResJ][idxResI][m][0][0], m);
				unsigned char r = mat[idxResJ][idxResI][m][0][0];
				unsigned char g = mat[idxResJ][idxResI][m][0][1];
				unsigned char b = mat[idxResJ][idxResI][m][0][2];
				//unsigned char r = mat[idxResJ][idxResI][m][0][0];
				if (suc != 1) {
					TIFFClose(tiff);
					return 0;
				}
			}
			TIFFClose(tiff);
			bear::dynamic_image_ptr image;
			try {
				image = bear::dynamic_image_ptr(w, h, cn, bear::data_type(2), depth / 8, &mat[idxResJ][idxResI][0][0][0], w*cn*depth / 8);
			}
			catch (bear::bear_exception e) {
				std::cout << e.what() << std::endl;
				return 0;
			}
			srcMat[idxResJ].push_back(image);

		}
	}
	width = width + 4;
	height = height + 4;
	bear::tensor<char, 3> dst;
	bear::dynamic_image_ptr dstImage;
	if (depth > 8) {
		dst = bear::tensor<char, 3>(height, width, cn*2);
	}
	else {
		dst = bear::tensor<char, 3>(height, width, cn);
	}
	try {
		dstImage = bear::dynamic_image_ptr(width, height, cn, bear::data_type(2), depth / 8, &dst[0][0][0], width*cn*depth / 8);
	}	
	catch (bear::bear_exception e) {
		std::cout << e.what() << std::endl;
		return 0;
	}

	//std::vector<std::vector<bear::dynamic_image_ptr>>* border = NULL;
	if (mode){
		vector<vector<bear::tensor<char, 3>>> mat(1, vector<bear::tensor< char, 3>>(param));
		vector<vector<bear::dynamic_image_ptr>> borderMat(1);
		for (int idxR = 0; idxR < param; idxR++)
		{
			string picFilePath = string(workDirPath) + "/_r_";
			picFilePath += to_string(idxR);
			picFilePath += "/channels/";
			picFilePath += fileName;

			TIFF* tiff = TIFFOpen(picFilePath.c_str(), "r");
			if (tiff == NULL) {
				return 0;
			}
			uint32 w, h;

			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);

			if (depth > 8) {
				mat[0][idxR] = bear::tensor<char, 3>(h, w, cn * 2);
			}
			else {
				mat[0][idxR] = bear::tensor<char, 3>(h, w, cn);
			}

			for (int m = 0; m < h; m++)
			{
				int suc = TIFFReadScanline(tiff, &mat[0][idxR][m][0][0], m);
				if (suc != 1) {
					return 0;
				}
			}
			bear::dynamic_image_ptr image;
			try {
				image = bear::dynamic_image_ptr(w, h, cn, bear::data_type(2), depth / 8, &mat[0][idxR][0][0][0], w*cn*depth / 8);
			}
			catch (bear::bear_exception e) {
				std::cout << e.what() << std::endl;
				return 0;
			}
			borderMat[0].push_back(image);
		}

		try {
			go_poisson_stiching(&dstImage, &srcMat, &borderMat, 2, 3, mode);
		//	image = bear::dynamic_image_ptr(w, h, cn, bear::data_type(2), depth / 8, &mat[0][idxR][0][0][0], width*cn*depth / 8);
		}
		catch (bear::bear_exception e) {
			std::cout << e.what() << std::endl;
			return 0;
		}
			

	}
	else {
		try {
			go_poisson_stiching(&dstImage, &srcMat, NULL, 2, 3, mode);
			//	image = bear::dynamic_image_ptr(w, h, cn, bear::data_type(2), depth / 8, &mat[0][idxR][0][0][0], width*cn*depth / 8);
		}
		catch (bear::bear_exception e) {
			std::cout << e.what() << std::endl;
			return 0;
		}
	}
	string resultPath = string(workDirPath) ;
	resultPath += "/final/";
	resultPath += fileName;
	TIFF* out = TIFFOpen(resultPath.c_str(), "w");





	uint16 compression = COMPRESSION_LZW;//
	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(out, TIFFTAG_COMPRESSION, compression);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, depth);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, cn);
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photoMetric);
	TIFFSetField(out, TIFFTAG_PLANARCONFIG, planarConfig);



	for (int m = 0; m < height; m++)
	{
		


		
		//if (m < height / 2) {
		//	suc = TIFFWriteScanline(out, buffer, m);
		//	if (suc != 1) {
		//		suc = 2;
		//	}
		//}
		unsigned char r = dstImage.data()[m*width*cn*depth / 8];
		unsigned char g = dstImage.data()[m*width*cn*depth / 8+1];
		unsigned char b = dstImage.data()[m*width*cn*depth / 8+2];

		int suc = TIFFWriteScanline(out, &dstImage.data()[m*width*cn*depth/8], m);
		if (suc != 1) {
			return 0;
		}

	}
	
	TIFFClose(out);



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
