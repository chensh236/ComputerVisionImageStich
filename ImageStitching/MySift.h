#ifdef _CH_  
#pragma package <opencv>  
#endif  

#ifndef _EiC  
#include <stdio.h>  

#include "stdlib.h"  
#include "string.h"   
#include "malloc.h"   
#include "math.h"   
#include <assert.h>  
#include <ctype.h>  
#include <time.h>  
//#include <cv.h>
//#include <cxcore.h>
#include <CImg.h>
//#include <highgui.h>
#include <vector>  
#include <iostream>
using namespace std;
#endif  

#ifdef _EiC  
#define WIN32  
#endif  

#define NUMSIZE 2  
#define GAUSSKERN 3.5  
#define PI 3.14159265358979323846  
#define null 0

//Sigma of base image -- See D.L.'s paper.  
#define INITSIGMA 0.5  
//Sigma of each octave -- See D.L.'s paper.  
#define SIGMA sqrt(3)//1.6//  

//Number of scales per octave.  See D.L.'s paper.  
#define SCALESPEROCTAVE 2  
#define MAXOCTAVES 4  

#define CONTRAST_THRESHOLD   0.02  
#define CURVATURE_THRESHOLD  10.0  
#define DOUBLE_BASE_IMAGE_SIZE 1  
#define peakRelThresh 0.8  
#define LEN 128 
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

//特征描述点，网格    
#define GridSpacing 4
using namespace cimg_library;

//Data structure for a float image.  
typedef struct ImageSt {        /*金字塔每一层*/
	float levelsigma;
	int levelsigmalength;
	float absolute_sigma;
	CImg<float> Level;       //CvMat是OPENCV的矩阵类，其元素可以是图像的象素值
} ImageLevels;

typedef struct ImageSt1 {      /*金字塔每一阶梯*/
	int row, col;          //Dimensions of image.   
	float subsample;
	ImageLevels Octave[5];
} ImageOctaves;

//keypoint数据结构，Lists of keypoints are linked by the "next" field.  
typedef struct KeypointSt {
	float row, col;     /* 反馈回原图像大小，特征点的位置 */
	float sx, sy;       /* 金字塔中特征点的位置 */
	int octave, level;  /* 金字塔中，特征点所在的阶梯、层次 */

	float scale, ori, mag;    /* 所在层的尺度sigma,主方向orientation (range [-PI,PI])，以及幅值 */
	float *descrip;           /* 特征描述字指针：128维或32维等 */
//	struct KeypointSt *next;  /* Pointer to next keypoint in list. */
} Keypoint;

class MySift {
public:
	MySift();
	~MySift();

	MySift(char* _filename, int _isColor);

	CImg<float> halfSizeImage(CImg<float> im);     //缩小图像：下采样  
	CImg<float> doubleSizeImage(CImg<float> im);   //扩大图像：最近临方法  
	CImg<float> doubleSizeImage2(CImg<float> im);  //扩大图像：线性插值  
	float getPixelBI(CImg<float> im, float col, float row);//双线性插值函数  
	void normalizeVec(float* vec, int dim);//向量归一化    
	CImg<float> GaussianKernel2D(float sigma);  //得到2维高斯核  
	void normalizeMat(CImg<float>& mat);        //矩阵归一化


    CImg<float> useFilter(CImg<float> img_in, vector<float> filterIn);
    vector<float> GaussianKernel1D(int dim, float sigma);
	//float* GaussianKernel1D(float sigma, int dim); //得到1维高斯核

//	//在具体像素处宽度方向进行高斯卷积
//	float ConvolveLocWidth(float* kernel, int dim, CImg<float> src, int x, int y);
//	//在整个图像宽度方向进行1D高斯卷积
//	void Convolve1DWidth(float* kern, int dim, CImg<float> src, CImg<float>& dst);
//	//在具体像素处高度方向进行高斯卷积
//	float ConvolveLocHeight(float* kernel, int dim, CImg<float> src, int x, int y);
//	//在整个图像高度方向进行1D高斯卷积
//	void Convolve1DHeight(float* kern, int dim, CImg<float> src, CImg<float>& dst);
	//用高斯函数模糊图像    
	int BlurImage(CImg<float> src, CImg<float>& dst, float sigma);


	//SIFT算法第一步：图像预处理  
	CImg<float> ScaleInitImage(CImg<float> im);                  //金字塔初始化  

	//SIFT算法第二步：建立高斯金字塔函数  
	void BuildGaussianOctaves(CImg<float> image);  //建立高斯金字塔

	//SIFT算法第三步：特征点位置检测，最后确定特征点的位置  
	void DetectKeypoint(int numoctaves, ImageOctaves *GaussianPyr);
	void DisplayKeypointLocation(CImg<float>& image, ImageOctaves *GaussianPyr);

	//SIFT算法第四步：计算高斯图像的梯度方向和幅值，计算各个特征点的主方向  
	void ComputeGrad_DirecandMag(int numoctaves, ImageOctaves *GaussianPyr);

	int FindClosestRotationBin(int binCount, float angle);  //进行方向直方图统计  
	void AverageWeakBins(float* bins, int binCount);       //对方向直方图滤波  
															//确定真正的主方向  
	bool InterpolateOrientation(float left, float middle, float right, float *degreeCorrection, float *peakValue);
	//确定各个特征点处的主方向函数  
	void AssignTheMainOrientation(int numoctaves, ImageOctaves *GaussianPyr, ImageOctaves *mag_pyr, ImageOctaves *grad_pyr);
	//显示主方向  
	void DisplayOrientation(CImg<float> image, ImageOctaves *GaussianPyr);

	//SIFT算法第五步：抽取各个特征点处的特征描述字  
	void ExtractFeatureDescriptors(int numoctaves, ImageOctaves *GaussianPyr);

	//为了显示图象金字塔，而作的图像水平、垂直拼接  
	CImg<float> MosaicHorizen(CImg<float> im1, CImg<float> im2);
	CImg<float> MosaicVertical(CImg<float> im1, CImg<float> im2);

	/* 以下为在源代码基础上添加部分 */
	void SiftMainProcess();
	vector<Keypoint> getFirstKeyDescriptors();  //获取第一个keyDescriptor结点

	void saveImgWithKeypoint(char* filename);
    float ImLevelsDog(int i, int j, int row, int col);
    float ImLevelsGussian(int i, int j, int row, int col);
    CImg<float> toGrayScale(CImg<float> img);
    void cvConvertScale(CImg<float> src, CImg<float> &dst, float scale, float shift);
    void cvConvertScaleAbs(CImg<float> src, CImg<float> &dst, float scale, float shift);
    void cvMinMaxLoc(CImg<float> img, float *min, float *max);
    void cvAddS(CImg<float> img, float input, CImg<float> &dst);
	void cvSub(CImg<float> im1, CImg<float> im2, CImg<float> &dst);

private:
	char* filename;
	int isColor;

	int numoctaves;

	ImageOctaves DOGoctaves[4];
	//DOG pyr，DOG算子计算简单，是尺度归一化的LoG算子的近似。  

	ImageOctaves *mag_thresh;
	ImageOctaves mag_pyr[4];
	ImageOctaves grad_pyr[4];

	int keypoint_count = 0;

	//定义特征点具体变量  
	vector<Keypoint> keypoints;      //用于临时存储特征点的位置等
	vector<Keypoint> keyDescriptors; //用于最后的确定特征点以及特征描述字

	const unsigned char color[3] = {255, 0, 0};
	//声明当前帧IplImage指针  
	CImg<float> src;
	CImg<float> image_kp;
	CImg<float> image_featDir;
	CImg<float> grey_im1;

	CImg<float> mosaic1;
	CImg<float> mosaic2;

	CImg<float> mosaicHorizen1;
	CImg<float> mosaicHorizen2;
	CImg<float> mosaicVertical1;

	CImg<float> image1Mat;
	CImg<float> tempMat;

	ImageOctaves Gaussianpyr[4];
};