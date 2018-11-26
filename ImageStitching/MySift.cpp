#include "MySift.h"
#include <float.h>

MySift::MySift() {
	keyDescriptors.clear();
	keyDescriptors.clear();
}

MySift::~MySift() {
}

MySift::MySift(char* _filename, int _isColor) {
	filename = _filename;
	isColor = _isColor;
}

//下采样原来的图像，返回缩小2倍尺寸的图像  
CImg<float> MySift::halfSizeImage(CImg<float> im) {
	int w = im.width() / 2;
	int h = im.height() / 2;
	CImg<float> imnew(w, h, 1, 1, 0);
    cimg_forXY(imnew, x, y){
        imnew(x, y) = im(x * 2, y * 2);
    }
	return imnew;
}

//上采样原来的图像，返回放大2倍尺寸的线性插值图像  
CImg<float> MySift::doubleSizeImage2(CImg<float> im) {
    int w = im.width() * 2;
    int h = im.height() * 2;
    CImg<float> imnew(w, h, 1, 1, 0);

    cimg_forXY(imnew, x, y){
        imnew(x, y) = im(x / 2, y / 2);
    }
    /*
    A B C
    E F G
    H I J
    pixels A C H J are pixels from original image
    pixels B E G I F are interpolated pixels
    */
    // interpolate pixels B and I
    int i, j;
    for ( j = 0; j < h; j += 2)
        for ( i = 1; i < w - 1; i += 2)
            imnew(i,j)=0.5*(im(i/2, j/2)+im(i/2+1, j/2));
    // interpolate pixels E and G
    for ( j = 1; j < h - 1; j += 2)
        for ( i = 0; i < w; i += 2)
            imnew(i,j)=0.5*(im(i/2, j/2)+im(i/2, j/2+1));
    // interpolate pixel F
    for ( j = 1; j < h - 1; j += 2)
        for ( i = 1; i < w - 1; i += 2)
            imnew(i,j)=0.25*(im(i/2, j/2)+im(i/2, j/2+1)+im(i/2+1, j/2)+im(i/2+1, j/2+1));
    return imnew;
}

//双线性插值，返回像素间的灰度值  
float MySift::getPixelBI(CImg<float> im, float col, float row) {
	int irow, icol;
	float rfrac, cfrac;
	float row1 = 0, row2 = 0;
	int width = im.width();
	int height = im.height(); 

	irow = (int)row;
	icol = (int)col;

	if (irow < 0 || irow >= height
		|| icol < 0 || icol >= width)
		return 0;
	if (row > height - 1)
		row = height - 1;
	if (col > width - 1)
		col = width - 1;
	rfrac = 1.0 - (row - (float)irow);
	cfrac = 1.0 - (col - (float)icol);
	if (cfrac < 1) {
		row1 = cfrac * im(icol, irow) + (1.0 - cfrac) * im(icol + 1, irow);
	}
	else {
		row1 = im(icol, irow);
	}
	if (rfrac < 1) {
		if (cfrac < 1) {
			row2 = cfrac * im(icol, irow + 1) + (1.0 - cfrac) * im(icol + 1, irow + 1);
		}
		else {
			row2 = im(icol, irow + 1);
		}
	}
	return rfrac * row1 + (1.0 - rfrac) * row2;
}

//向量归一化
void MySift::normalizeVec(float* vec, int dim)
{
	unsigned int i;
	float sum = 0;
	for (i = 0; i < dim; i++)
		sum += vec[i];
	for (i = 0; i < dim; i++)
		vec[i] /= sum;
}

//得到向量的欧式长度，2-范数  
float GetVectorNorm(float* vec, int dim)
{
	float sum = 0.0;
	for (unsigned int i = 0; i<dim; i++)
		sum += vec[i] * vec[i];
	return sqrt(sum);
}

//产生2D高斯核矩阵  
CImg<float> MySift::GaussianKernel2D(float sigma)
{
	// int dim = (int) max(3.0f, GAUSSKERN * sigma);  
	int dim = (int)max(3.0f, 2.0 * GAUSSKERN *sigma + 1.0f);
	// make dim odd  
	if (dim % 2 == 0)
		dim++;
	//printf("GaussianKernel(): Creating %dx%d matrix for sigma=%.3f gaussian/n", dim, dim, sigma);  
	CImg<float> mat(dim, dim, 1, 1, 0);
	float s2 = sigma * sigma;
	int c = dim / 2;
	//printf("%d %d/n", mat.size(), mat[0].size());  
	//!
    float m = 1.0 / (sqrt(2.0 * PI) * sigma);
	for (int i = 0; i < (dim + 1) / 2; i++)
	{
		for (int j = 0; j < (dim + 1) / 2; j++)
		{
			//printf("%d %d %d/n", c, i, j);  
			float v = m * exp(-(1.0*i*i + 1.0*j*j) / (2.0 * s2));
			mat(c + i, c + j) = v;
			mat(c - i, c + j) = v;
			mat(c + i, c - j) = v;
			mat(c - i, c - j) = v;
		}
	}
	// normalizeMat(mat);  
	return mat;
}

vector<float> MySift::GaussianKernel1D(int dim, float sigma)
{
    vector<float> kern;
    for(int i = 0; i < dim; ++ i){
        kern.push_back(0);
    }
    float s2 = sigma * sigma;
    int c = dim / 2;
    float m= 1.0/(sqrt(2.0 * PI) * sigma);
    double v;
    for (int  i = 0; i < (dim + 1) / 2; i++)
    {
        v = m * exp(-(1.0*i*i)/(2.0 * s2)) ;
        kern[c + i] = v;
        kern[c - i] = v;
    }
    //   normalizeVec(kern, dim);
    // for ( i = 0; i < dim; i++)
    //  printf("%f  ", kern[i]);
    //  printf("/n");
    return kern;
}

CImg<float> MySift::useFilter(CImg<float> img_in, vector<float> filterIn)
{
	int dim = filterIn.size();
	int cen = dim / 2;
	CImg<float> Xresult(img_in.width(), img_in.height(), 1, 1, 0);
	CImg<float> result = Xresult;

	// x direction
	cimg_forXY(Xresult, x, y){
		for(int i = 0; i < dim; ++i){
			int col = x + i - cen;
			if(col < 0) col = 0;
			if(col >= img_in.width()) col = img_in.width() - 1;
			Xresult(x, y) += filterIn[i] * img_in(col, y);
		}
		if(Xresult(x, y) > 1) Xresult(x, y) = 1;
	}

	cimg_forXY(result, x, y){
		for(int i = 0; i < dim; i++){
			int row = y + i - cen;
			if(row < 0) row = 0;
			if(row >= img_in.height()) row = img_in.height() - 1;
			result(x, y) += filterIn[i] * Xresult(x, row);
		}
		if(result(x, y) > 1) result(x, y) = 1;
	}
	return result;
}


//卷积模糊图像  
int MySift::BlurImage(CImg<float> src, CImg<float>& dst, float sigma)
{
	vector<float> convkernel;
	int dim = (int)max(3.0f, 2.0 * GAUSSKERN * sigma + 1.0f);
	// make dim odd  
	if (dim % 2 == 0)
		dim++;
	convkernel = GaussianKernel1D(dim, sigma);

	dst = useFilter(src, convkernel);
	return dim;
}

CImg<float> MySift::ScaleInitImage(CImg<float> im)
{
	float sigma, preblur_sigma;
	
	CImg<float>tempMat;

	//首先对图像进行平滑滤波，抑制噪声  
	CImg<float>imMat(im.width(), im.height(), 1, 1, 0);
    CImg<float>dst(im.width(), im.height(), 1, 1, 0);
	BlurImage(im, imMat, INITSIGMA);

	//针对两种情况分别进行处理：初始化放大原始图像或者在原图像基础上进行后续操作  
	//建立金字塔的最底层  
	if (DOUBLE_BASE_IMAGE_SIZE)
	{
		tempMat = doubleSizeImage2(imMat);//对扩大两倍的图像进行二次采样，采样率为0.5，采用线性插值
		preblur_sigma = 1.0;//sqrt(2 - 4*INITSIGMA*INITSIGMA);  
		BlurImage(tempMat, dst, preblur_sigma);
		// The initial blurring for the first image of the first octave of the pyramid.  
		sigma = sqrt((4 * INITSIGMA*INITSIGMA) + preblur_sigma * preblur_sigma);
		//  sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA * 4);  
		//printf("Init Sigma: %f/n", sigma);  
		BlurImage(dst, tempMat, sigma);       //得到金字塔的最底层-放大2倍的图像  
		//cvReleaseMat(&dst);
		return tempMat;
	}
	else
	{
		//sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA);  
		preblur_sigma = 1.0;//sqrt(2 - 4*INITSIGMA*INITSIGMA);  
		sigma = sqrt((4 * INITSIGMA*INITSIGMA) + preblur_sigma * preblur_sigma);
		//printf("Init Sigma: %f/n", sigma);  
		BlurImage(imMat, dst, sigma);        //得到金字塔的最底层：原始图像大小  
		return dst;
	}
}

void MySift::cvSub(CImg<float> im1, CImg<float> im2, CImg<float> &dst){
    cimg_forXY(dst, x, y){
        dst(x, y) = im1(x, y) - im2(x, y);
    }
}

//SIFT算法第二步
void MySift::BuildGaussianOctaves(CImg<float> image)
{
	CImg<float>tempMat;
	CImg<float>dst;
	CImg<float>temp;

	int i, j;
	float k = pow(2, 1.0 / ((float)SCALESPEROCTAVE));  //方差倍数  
	float preblur_sigma, initial_sigma, sigma1, sigma2, sigma, absolute_sigma, sigma_f;
	//计算金字塔的阶梯数目  
	int dim = min(image.height(), image.width());
	int numoctaves = (int)(log((float)dim) / log(2.0)) - 2;    //金字塔阶数  
																//限定金字塔的阶梯数  
	numoctaves = min(numoctaves, MAXOCTAVES);
	//为高斯金塔和DOG金字塔分配内存

	printf("BuildGaussianOctaves(): Base image dimension is %dx%d\n", (int)(0.5*(image.width())), (int)(0.5*(image.height())));
	printf("BuildGaussianOctaves(): Building %d octaves\n", numoctaves);

	// start with initial source image
	// temp是最底层
	tempMat = image;
	// preblur_sigma = 1.0;//sqrt(2 - 4*INITSIGMA*INITSIGMA);  
	initial_sigma = sqrt(2.0);//sqrt( (4*INITSIGMA*INITSIGMA) + preblur_sigma * preblur_sigma );  
							  //   initial_sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA * 4);  

							  //在每一阶金字塔图像中建立不同的尺度图像  
	for (i = 0; i < numoctaves; i++)
	{
		//首先建立金字塔每一阶梯的最底层，其中0阶梯的最底层已经建立好  
		printf("Building octave %d of dimesion (%d, %d)\n", i, tempMat.width(), tempMat.height());
		//存储各个阶梯的最底层  
		(Gaussianpyr[i].Octave)[0].Level = tempMat;

		Gaussianpyr[i].col = tempMat.width();
		Gaussianpyr[i].row = tempMat.height();
		DOGoctaves[i].col = tempMat.width();
		DOGoctaves[i].row = tempMat.height();
		if (DOUBLE_BASE_IMAGE_SIZE)
			Gaussianpyr[i].subsample = pow(2.0, i)*0.5;
		else
			Gaussianpyr[i].subsample = pow(2.0, i);

		if (i == 0)
		{
			(Gaussianpyr[0].Octave)[0].levelsigma = initial_sigma;
			(Gaussianpyr[0].Octave)[0].absolute_sigma = initial_sigma;
			printf("0 scale and blur sigma : %f \n", (Gaussianpyr[0].subsample) * ((Gaussianpyr[0].Octave)[0].absolute_sigma));
		}
		else
		{
			(Gaussianpyr[i].Octave)[0].levelsigma = (Gaussianpyr[i - 1].Octave)[SCALESPEROCTAVE].levelsigma;
			(Gaussianpyr[i].Octave)[0].absolute_sigma = (Gaussianpyr[i - 1].Octave)[SCALESPEROCTAVE].absolute_sigma;
			printf("0 scale and blur sigma : %f \n", ((Gaussianpyr[i].Octave)[0].absolute_sigma));
		}
		sigma = initial_sigma;
		//建立本阶梯其他层的图像  
		for (j = 1; j < SCALESPEROCTAVE + 3; j++)
		{
			dst = CImg<float>(tempMat.width(), tempMat.height(), 1, 1, 0);//用于存储高斯层  
			temp = CImg<float>(tempMat.width(), tempMat.height(), 1, 1, 0);//用于存储DOG层  
																	   // 2 passes of 1D on original  
																	   //   if(i!=0)  
																	   //   {  
																	   //       sigma1 = pow(k, j - 1) * ((Gaussianpyr[i-1].Octave)[j-1].levelsigma);  
																	   //          sigma2 = pow(k, j) * ((Gaussianpyr[i].Octave)[j-1].levelsigma);  
																	   //       sigma = sqrt(sigma2*sigma2 - sigma1*sigma1);  
			sigma_f = sqrt(k*k - 1)*sigma;
			//   }  
			//   else  
			//   {  
			//       sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA * 4)*pow(k,j);  
			//   }    
			sigma = k*sigma;
			absolute_sigma = sigma * (Gaussianpyr[i].subsample);
			printf("%d scale and Blur sigma: %f  \n", j, absolute_sigma);

			(Gaussianpyr[i].Octave)[j].levelsigma = sigma;
			(Gaussianpyr[i].Octave)[j].absolute_sigma = absolute_sigma;
			//产生高斯层  
			int length = BlurImage((Gaussianpyr[i].Octave)[j - 1].Level, dst, sigma_f);//相应尺度  
			(Gaussianpyr[i].Octave)[j].levelsigmalength = length;
			(Gaussianpyr[i].Octave)[j].Level = dst;
			//产生DOG层  
			cvSub(((Gaussianpyr[i].Octave)[j]).Level, ((Gaussianpyr[i].Octave)[j - 1]).Level, temp);
			//         cvAbsDiff( ((Gaussianpyr[i].Octave)[j]).Level, ((Gaussianpyr[i].Octave)[j-1]).Level, temp );  
			((DOGoctaves[i].Octave)[j - 1]).Level = temp;
		}
		// halve the image size for next iteration  
		tempMat = halfSizeImage(((Gaussianpyr[i].Octave)[SCALESPEROCTAVE].Level));
	}
}

float MySift::ImLevelsDog(int i, int j, int row, int col){
    return (DOGoctaves[i].Octave)[j].Level(col, row);
}

float MySift::ImLevelsGussian(int i, int j, int row, int col){
    return (Gaussianpyr[i].Octave)[j].Level(col, row);
}


//SIFT算法第三步，特征点位置检测，  
void MySift::DetectKeypoint(int numoctaves, ImageOctaves *GaussianPyr) {
	//计算用于DOG极值点检测的主曲率比的阈值  
	float curvature_threshold;
	curvature_threshold = ((CURVATURE_THRESHOLD + 1)*(CURVATURE_THRESHOLD + 1)) / CURVATURE_THRESHOLD;
	// 遍历所有的金字塔
    for (int i = 0; i < numoctaves; i++) {
		for (int j = 1; j < SCALESPEROCTAVE + 1; j++) {    //取中间的scaleperoctave个层
														   //在图像的有效区域内寻找具有显著性特征的局部最大值  
														   //float sigma=(GaussianPyr[i].Octave)[j].levelsigma;  
														   //int dim = (int) (max(3.0f, 2.0*GAUSSKERN *sigma + 1.0f)*0.5);  
			int dim = (int)(0.5*((GaussianPyr[i].Octave)[j].levelsigmalength) + 0.5);
			// 寻找领域
			for (int m = dim; m < ((DOGoctaves[i].row) - dim); m++)
				for (int n = dim; n < ((DOGoctaves[i].col) - dim); n++) {
					if (fabs(ImLevelsDog(i, j, m, n)) >= CONTRAST_THRESHOLD) {

						if (ImLevelsDog(i, j, m, n) != 0.0) {    //1、首先是非零  
							// 对应的该值
							float inf_val = ImLevelsDog(i, j, m, n);
							// 极小值
							// 上一层的9个点
							if (((inf_val <= ImLevelsDog(i, j - 1, m - 1, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m + 1, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m - 1, n)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m, n)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m + 1, n)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m - 1, n + 1)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m, n + 1)) &&
								(inf_val <= ImLevelsDog(i, j - 1, m + 1, n + 1)) &&    //底层的小尺度9  
								// 该层的8个点
								(inf_val <= ImLevelsDog(i, j, m - 1, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j, m, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j, m + 1, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j, m - 1, n)) &&
								(inf_val <= ImLevelsDog(i, j, m + 1, n)) &&
								(inf_val <= ImLevelsDog(i, j, m - 1, n + 1)) &&
								(inf_val <= ImLevelsDog(i, j, m, n + 1)) &&
								(inf_val <= ImLevelsDog(i, j, m + 1, n + 1)) &&     //当前层8  
								// 下一层的9个点
								(inf_val <= ImLevelsDog(i, j + 1, m - 1, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m + 1, n - 1)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m - 1, n)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m, n)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m + 1, n)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m - 1, n + 1)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m, n + 1)) &&
								(inf_val <= ImLevelsDog(i, j + 1, m + 1, n + 1))     //下一层大尺度9          
								) ||
								// 极大值
								// 上一层的9个点
								((inf_val >= ImLevelsDog(i, j - 1, m - 1, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m + 1, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m - 1, n)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m, n)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m + 1, n)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m - 1, n + 1)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m, n + 1)) &&
									(inf_val >= ImLevelsDog(i, j - 1, m + 1, n + 1)) &&
									// 该层的8个点
									(inf_val >= ImLevelsDog(i, j, m - 1, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j, m, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j, m + 1, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j, m - 1, n)) &&
									(inf_val >= ImLevelsDog(i, j, m + 1, n)) &&
									(inf_val >= ImLevelsDog(i, j, m - 1, n + 1)) &&
									(inf_val >= ImLevelsDog(i, j, m, n + 1)) &&
									(inf_val >= ImLevelsDog(i, j, m + 1, n + 1)) &&
									// 下一层的9个点
									(inf_val >= ImLevelsDog(i, j + 1, m - 1, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m + 1, n - 1)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m - 1, n)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m, n)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m + 1, n)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m - 1, n + 1)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m, n + 1)) &&
									(inf_val >= ImLevelsDog(i, j + 1, m + 1, n + 1))
									)) {           //2、满足26个中极值点  

												   //此处可存储  
												   //然后必须具有明显的显著性，即必须大于CONTRAST_THRESHOLD=0.02  
												   //去除相应较弱的点
								if (fabs(ImLevelsDog(i, j, m, n)) >= CONTRAST_THRESHOLD) {
									//最后显著处的特征点必须具有足够的曲率比，CURVATURE_THRESHOLD=10.0，首先计算Hessian矩阵  
									// Compute the entries of the Hessian matrix at the extrema location.  
									/*
									1   0   -1
									0   0   0
									-1   0   1         *0.25
									*/
									// Compute the trace and the determinant of the Hessian.  
									//Tr_H = Dxx + Dyy;  
									//Det_H = Dxx*Dyy - Dxy^2;  
									float Dxx, Dyy, Dxy, Tr_H, Det_H, curvature_ratio;
									Dxx = ImLevelsDog(i, j, m, n - 1) + ImLevelsDog(i, j, m, n + 1) - 2.0*ImLevelsDog(i, j, m, n);
									Dyy = ImLevelsDog(i, j, m - 1, n) + ImLevelsDog(i, j, m + 1, n) - 2.0*ImLevelsDog(i, j, m, n);
									Dxy = ImLevelsDog(i, j, m - 1, n - 1) + ImLevelsDog(i, j, m + 1, n + 1) - ImLevelsDog(i, j, m + 1, n - 1) - ImLevelsDog(i, j, m - 1, n + 1);
									Tr_H = Dxx + Dyy;
									Det_H = Dxx*Dyy - Dxy*Dxy;
									// Compute the ratio of the principal curvatures.  
									// PPT 36
									curvature_ratio = (1.0*Tr_H*Tr_H) / Det_H;
									if ((Det_H >= 0.0) && (curvature_ratio <= curvature_threshold)) {    //最后得到最具有显著性特征的特征点  

																										 //将其存储起来，以计算后面的特征描述字
										Keypoint k;
										k.row = m*(GaussianPyr[i].subsample);
										k.col = n*(GaussianPyr[i].subsample);
										k.sy = m;    //行  
										k.sx = n;    //列  
										k.octave = i;
										k.level = j;
										// 所在层的尺度
										k.scale = (GaussianPyr[i].Octave)[j].absolute_sigma;
										keypoints.push_back(k);
									}//if >curvature_thresh  
								}//if >contrast  
							}//if inf value  极值
						}//if non zero  
					}//if >contrast  
				}  //for concrete image level col  
		}//for levels  
	}//for octaves
}

//在图像中，显示SIFT特征点的位置  
void MySift::DisplayKeypointLocation(CImg<float>& image, ImageOctaves *GaussianPyr) {
	for(Keypoint p : keypoints) {
        image.draw_circle((int)(p.col), (int)(p.row), 2, color);
	}
}

// Compute the gradient direction and magnitude of the gaussian pyramid images  
void MySift::ComputeGrad_DirecandMag(int numoctaves, ImageOctaves *GaussianPyr)
{
	// ImageOctaves *mag_thresh ;
	// float sigma=( (GaussianPyr[0].Octave)[SCALESPEROCTAVE+2].absolute_sigma ) / GaussianPyr[0].subsample;  
	// int dim = (int) (max(3.0f, 2 * GAUSSKERN *sigma + 1.0f)*0.5+0.5);  
	for (int i = 0; i<numoctaves; i++)
	{
		for (int j = 1; j<SCALESPEROCTAVE + 1; j++)//取中间的scaleperoctave个层  
		{
			CImg<float>Mag(GaussianPyr[i].col, GaussianPyr[i].row, 1, 1, 0);
			CImg<float>Ori(GaussianPyr[i].col, GaussianPyr[i].row, 1, 1, 0);
			CImg<float>tempMat1(GaussianPyr[i].col, GaussianPyr[i].row, 1, 1, 0);
			CImg<float>tempMat2(GaussianPyr[i].col, GaussianPyr[i].row, 1, 1, 0);
			
			for (int m = 1; m<(GaussianPyr[i].row - 1); m++)
				for (int n = 1; n<(GaussianPyr[i].col - 1); n++)
				{
					//计算幅值
					// PPT 38 差分
					tempMat1(n, m) = 0.5*(ImLevelsGussian(i, j, m, n + 1) - ImLevelsGussian(i, j, m, n - 1));  //dx  
					tempMat2(n, m) = 0.5*(ImLevelsGussian(i, j, m + 1, n) - ImLevelsGussian(i, j, m - 1, n));  //dy  
					Mag(n, m) = sqrt(tempMat1(n, m)*tempMat1(n, m) + tempMat2(n, m)*tempMat2(n, m));  //mag
																									  //计算方向
					Ori(n, m) = atan(tempMat2(n, m) / tempMat1(n, m));
					if (Ori(n, m) == PI/**/)
						Ori(n, m) = -PI/**/;
				}
			((mag_pyr[i].Octave)[j - 1]).Level = Mag;
			((grad_pyr[i].Octave)[j - 1]).Level = Ori;
		}//for levels  
	}//for octaves  
}

//SIFT算法第四步：计算各个特征点的主方向，确定主方向  
void MySift::AssignTheMainOrientation(int numoctaves, ImageOctaves *GaussianPyr, ImageOctaves *mag_pyr, ImageOctaves *grad_pyr)
{
	// Set up the histogram bin centers for a 36 bin histogram.  
	int num_bins = 36;
	float hist_step = 2.0 * PI / num_bins;
	float hist_orient[36];
	// 所有的梯度
	for (int i = 0; i<36; i++)
		hist_orient[i] = -PI + i*hist_step;
	float sigma1 = (((GaussianPyr[0].Octave)[SCALESPEROCTAVE].absolute_sigma)) / (GaussianPyr[0].subsample);//SCALESPEROCTAVE+2  
	int zero_pad = (int)(max(3.0f, 2 * GAUSSKERN *sigma1 + 1.0f)*0.5 + 0.5);
	//Assign orientations to the keypoints.  

	for(Keypoint p : keypoints)
	{
		int i = p.octave;
		int j = p.level;
		int m = p.sy;   //行  
		int n = p.sx;   //列  
		if ((m >= zero_pad) && (m<GaussianPyr[i].row - zero_pad) &&
			(n >= zero_pad) && (n<GaussianPyr[i].col - zero_pad))
		{
			float sigma = (((GaussianPyr[i].Octave)[j].absolute_sigma)) / (GaussianPyr[i].subsample);
			//产生二维高斯模板  
			CImg<float> mat = GaussianKernel2D(sigma);
			int dim = (int)(0.5 * (mat.height()));
			//分配用于存储Patch幅值和方向的空间  

			//声明方向直方图变量  
			float* orienthist = (float *)malloc(36 * sizeof(float));
			for (int sw = 0; sw < 36; ++sw)
			{
				orienthist[sw] = 0.0;
			}

			//在特征点的周围统计梯度方向  
			for (int x = m - dim, mm = 0; x <= (m + dim); x++, mm++)
				for (int y = n - dim, nn = 0; y <= (n + dim); y++, nn++)
				{
					//计算特征点处的幅值  
					float dx = 0.5*(ImLevelsGussian(i, j, x, y + 1) - ImLevelsGussian(i, j, x, y - 1));  //dx
					float dy = 0.5*(ImLevelsGussian(i, j, x + 1, y) - ImLevelsGussian(i, j, x - 1, y));  //dy
					float mag = sqrt(dx*dx + dy*dy);  //mag  
													   //计算方向  
					float Ori = atan(1.0*dy / dx);
					int binIdx = FindClosestRotationBin(36, Ori);                   //得到离现有方向最近的直方块  
					orienthist[binIdx] = orienthist[binIdx] + 1.0* mag * mat(nn, mm);//利用高斯加权累加进直方图相应的块  
				}
				// 已经得到36个方向的块
			// Find peaks in the orientation histogram using nonmax suppression.  
			AverageWeakBins(orienthist, 36);
			// find the maximum peak in gradient orientation  
			float maxGrad = 0.0;
			// pos
			int maxBin = 0;
			for (int b = 0; b < 36; ++b)
			{
				if (orienthist[b] > maxGrad)
				{
					maxGrad = orienthist[b];
					maxBin = b;
				}
			}
			// First determine the real interpolated peak high at the maximum bin  
			// position, which is guaranteed to be an absolute peak.  
			float maxPeakValue = 0.0;
			float maxDegreeCorrection = 0.0;
			if ((InterpolateOrientation(orienthist[maxBin == 0 ? (36 - 1) : (maxBin - 1)],
				orienthist[maxBin], orienthist[(maxBin + 1) % 36],
				&maxDegreeCorrection, &maxPeakValue)) == false)
				printf("BUG: Parabola fitting broken");

			// Now that we know the maximum peak value, we can find other keypoint  
			// orientations, which have to fulfill two criterias:  
			//  
			//  1. They must be a local peak themselves. Else we might add a very  
			//     similar keypoint orientation twice (imagine for example the  
			//     values: 0.4 1.0 0.8, if 1.0 is maximum peak, 0.8 is still added  
			//     with the default threshhold, but the maximum peak orientation  
			//     was already added).  
			//  2. They must have at least peakRelThresh times the maximum peak  
			//     value.  
			bool binIsKeypoint[36];
			for (int b = 0; b < 36; ++b)
			{
				binIsKeypoint[b] = false;
				// The maximum peak of course is  
				if (b == maxBin)
				{
					binIsKeypoint[b] = true;
					continue;
				}
				// Local peaks are, too, in case they fulfill the threshhold  
				if (orienthist[b] < (peakRelThresh * maxPeakValue))
					continue;
				int leftI = (b == 0) ? (36 - 1) : (b - 1);
				int rightI = (b + 1) % 36;
				if (orienthist[b] <= orienthist[leftI] || orienthist[b] <= orienthist[rightI])
					continue; // no local peak  
				binIsKeypoint[b] = true;
			}
			// find other possible locations  
			float oneBinRad = (2.0 * PI) / 36;
			for (int b = 0; b < 36; ++b)
			{
				if (binIsKeypoint[b] == false)
					continue;
				int bLeft = (b == 0) ? (36 - 1) : (b - 1);
				int bRight = (b + 1) % 36;
				// Get an interpolated peak direction and value guess.  
				float peakValue;
				float degreeCorrection;

				float maxPeakValue, maxDegreeCorrection;
				if (InterpolateOrientation(orienthist[maxBin == 0 ? (36 - 1) : (maxBin - 1)],
					orienthist[maxBin], orienthist[(maxBin + 1) % 36],
					&degreeCorrection, &peakValue) == false)
				{
					printf("BUG: Parabola fitting broken");
				}

				float degree = (b + degreeCorrection) * oneBinRad - PI;
				if (degree < -PI)
					degree += 2.0 * PI;
				else if (degree > PI)
					degree -= 2.0 * PI;
				//存储方向，可以直接利用检测到的链表进行该步主方向的指定;  
				//分配内存重新存储特征点  
				Keypoint k;
				k.descrip = (float*)malloc(LEN * sizeof(float));
				k.row = p.row;
				k.col = p.col;
				k.sy = p.sy;    //行  
				k.sx = p.sx;    //列  
				k.octave = p.octave;
				k.level = p.level;
				k.scale = p.scale;
				k.ori = degree;
				k.mag = peakValue;
				keyDescriptors.push_back(k);
			}//for  
			free(orienthist);
		}
	}
    cout<<"points:"<<keypoints.size()<<endl;
}

//寻找与方向直方图最近的柱，确定其index   
int MySift::FindClosestRotationBin(int binCount, float angle)
{
	angle += PI/**/;
	angle /= 2.0 * PI/**/;
	// calculate the aligned bin  
	angle *= binCount;
	int idx = (int)angle;
	if (idx == binCount)
		idx = 0;
	return (idx);
}

// Average the content of the direction bins.
// 中值滤波
void MySift::AverageWeakBins(float* hist, int binCount)
{
	// TODO: make some tests what number of passes is the best. (its clear  
	// one is not enough, as we may have something like  
	// ( 0.4, 0.4, 0.3, 0.4, 0.4 ))  
	for (int sn = 0; sn < 2; ++sn)
	{
		float firstE = hist[0];
		float last = hist[binCount - 1];
		for (int sw = 0; sw < binCount; ++sw)
		{
			float cur = hist[sw];
			float next = (sw == (binCount - 1)) ? firstE : hist[(sw + 1) % binCount];
			hist[sw] = (last + cur + next) / 3.0;
			last = cur;
		}
	}
}

// Fit a parabol to the three points (-1.0 ; left), (0.0 ; middle) and  
// (1.0 ; right).  
// Formulas:  
// f(x) = a (x - c)^2 + b  
// c is the peak offset (where f'(x) is zero), b is the peak value.  
// In case there is an error false is returned, otherwise a correction  
// value between [-1 ; 1] is returned in 'degreeCorrection', where -1  
// means the peak is located completely at the left vector, and -0.5 just  
// in the middle between left and middle and > 0 to the right side. In  
// 'peakValue' the maximum estimated peak value is stored.  
// 二次曲线拟合
bool MySift::InterpolateOrientation(float left, float middle, float right, float *degreeCorrection, float *peakValue)
{
	float a = ((left + right) - 2.0 * middle) / 2.0;   //抛物线捏合系数a  
														// degreeCorrection = peakValue = float.NaN;  

														// Not a parabol  
	if (a == 0.0)
		return false;
	float c = (((left - middle) / a) - 1.0) / 2.0;
	float b = middle - c * c * a;
	if (c < -0.5 || c > 0.5)
		return false;
	*degreeCorrection = c;
	*peakValue = b;
	return true;
}

//显示特征点处的主方向  
void MySift::DisplayOrientation(CImg<float> image, ImageOctaves *GaussianPyr)
{
	for(Keypoint p : keyDescriptors) // 没到表尾
	{

		float scale = (GaussianPyr[p.octave].Octave)[p.level].absolute_sigma;
		float autoscale = 3.0;
		float uu = autoscale*scale*cos(p.ori);
		float vv = autoscale*scale*sin(p.ori);
		float x = (p.col) + uu;
		float y = (p.row) + vv;
        printf("%f %f %f %f\n",x, y, (p.col), (p.row));
        image.draw_line((int)(p.col), (int)(p.row), (int)x, (int)y, color);
		
//		// Arrow head parameters
		float alpha = 0.33; // Size of arrow head relative to the length of the vector
		float beta = 0.33;  // Width of the base of the arrow head relative to the length

		float xx0 = (p.col) + uu - alpha*(uu + beta*vv);
		float yy0 = (p.row) + vv - alpha*(vv - beta*uu);
		float xx1 = (p.col) + uu - alpha*(uu - beta*vv);
		float yy1 = (p.row) + vv - alpha*(vv + beta*uu);
        image.draw_line((int)xx0, (int)yy0, (int)x, (int)y, color);
        image.draw_line((int)xx1, (int)yy1, (int)x, (int)y, color);
	}
    image.display();
}

//SIFT算法第五步：抽取各个特征点处的特征描述字
void MySift::ExtractFeatureDescriptors(int numoctaves, ImageOctaves *GaussianPyr) {
	// The orientation histograms have 8 bins
	// 8个方向的角度值
	float orient_bin_spacing = PI / 4;
	float orient_angles[8] = { -PI, (float)-PI + orient_bin_spacing, -PI*0.5, -orient_bin_spacing,
		0.0, orient_bin_spacing, PI*0.5,  (float)PI + orient_bin_spacing };
	//产生描述字中心各点坐标
	//
	float *feat_grid = (float *)malloc(2 * 16 * sizeof(float));
	/*-6, -2   -------- 6,6
	 *
	 *
	 * */
	for (int i = 0; i < GridSpacing; i++) {
		for (int j = 0; j < 2 * GridSpacing; ++j, ++j) {
			feat_grid[i * 2 * GridSpacing + j] = -6.0 + i*GridSpacing;
			feat_grid[i * 2 * GridSpacing + j + 1] = -6.0 + 0.5*j*GridSpacing;
		}
	}

	/*
	 * -7.5 - 7.5 步长1
	 * */

	//产生网格  16个大网格，里面共有16个小网格
	float *feat_samples = (float *)malloc(2 * 256 * sizeof(float));
	for (int i = 0; i < 4 * GridSpacing; i++) {
		for (int j = 0; j < 8 * GridSpacing; j += 2) {
			feat_samples[i * 8 * GridSpacing + j] = -(2 * GridSpacing - 0.5) + i;
			feat_samples[i * 8 * GridSpacing + j + 1] = -(2 * GridSpacing - 0.5) + 0.5*j;
		}
	}
	for (int i = 0; i < 4 * GridSpacing; i++) {
		for (int j = 0; j < 8 * GridSpacing; j += 2) {
		}
	}


	float feat_window = 2 * GridSpacing;
	cout<<"descripsize:"<<keyDescriptors.size()<<endl;
	for(Keypoint p : keyDescriptors)
	{
		//float scale = (GaussianPyr[p.octave].Octave)[p.level].absolute_sigma;
		float sine = sin(p.ori);
		float cosine = cos(p.ori);

		//计算中心点坐标旋转之后的位置  
		float *featcenter = (float *)malloc(2 * 16 * sizeof(float));
		for (int i = 0; i < GridSpacing; i++) {
			for (int j = 0; j < 2 * GridSpacing; j += 2) {
				float x = feat_grid[i * 2 * GridSpacing + j];
				float y = feat_grid[i * 2 * GridSpacing + j + 1];
				featcenter[i * 2 * GridSpacing + j] = ((cosine * x + sine * y) + p.sx);
				featcenter[i * 2 * GridSpacing + j + 1] = ((-sine * x + cosine * y) + p.sy);
			}
		}

		// 網格中心點旋轉後的位置
		// calculate sample window coordinates (rotated along keypoint)  
		float *feat = (float *)malloc(2 * 256 * sizeof(float));
		for (int i = 0; i < 64 * GridSpacing; i++, i++) {
			float x = feat_samples[i];
			float y = feat_samples[i + 1];
			feat[i] = ((cosine * x + sine * y) + p.sx);
			feat[i + 1] = ((-sine * x + cosine * y) + p.sy);
		}

		// 初始化特征描述子
		//Initialize the feature descriptor.  
		float *feat_desc = (float *)malloc(LEN * sizeof(float));
		for (int i = 0; i < LEN; i++) {
			feat_desc[i] = 0.0;
			// printf("%f  ",feat_desc[i]);    
		}
		//printf("/n");  
		for (int i = 0; i < 512; ++i, ++i) {
			float x_sample = feat[i];
			float y_sample = feat[i + 1];

			// Interpolate the gradient at the sample position
			/*
                0   12   0
                21  22   23
                0   32   0   具体插值策略如图示
			*/
			float sample12 = getPixelBI(((GaussianPyr[p.octave].Octave)[p.level]).Level, x_sample, y_sample - 1);
			float sample21 = getPixelBI(((GaussianPyr[p.octave].Octave)[p.level]).Level, x_sample - 1, y_sample);
			//float sample22 = getPixelBI(((GaussianPyr[p.octave].Octave)[p.level]).Level, x_sample, y_sample);
			float sample23 = getPixelBI(((GaussianPyr[p.octave].Octave)[p.level]).Level, x_sample + 1, y_sample);
			float sample32 = getPixelBI(((GaussianPyr[p.octave].Octave)[p.level]).Level, x_sample, y_sample + 1);
			//float diff_x = 0.5*(sample23 - sample21);
			//float diff_y = 0.5*(sample32 - sample12);
			float diff_x = sample23 - sample21;
			float diff_y = sample32 - sample12;
			float mag_sample = sqrt(diff_x*diff_x + diff_y*diff_y);
			float grad_sample = atan(diff_y / diff_x);
			// 代替？
//			mag_sample = sqrt(diff_x*diff_x + diff_y*diff_y);
//			grad_sample = atan(diff_y / diff_x);
//			if (diff_x < 0)
//				grad_sample += CV_PI;
//			if (grad_sample >= CV_PI)   grad_sample -= 2*CV_PI;
			if (grad_sample == PI/**/)
				grad_sample = -PI/**/;
			//需要改進調優的地方！！！！！！！！！！
			// 計算採樣點對於4*4個種子點的權重，其實只有鄰近的種子點會有權重
			// 這類似 hog算子的 block 歸一化。
			// float[128]，對每個種子點內的8方向權值是一樣的。
			// Compute the weighting for the x and y dimensions.
			float *x_wght = (float *)malloc(GridSpacing * GridSpacing * sizeof(float));
			float *y_wght = (float *)malloc(GridSpacing * GridSpacing * sizeof(float));
			float *pos_wght = (float *)malloc(8 * GridSpacing * GridSpacing * sizeof(float));
			// 一共有16个点
			for (int m = 0; m < 32; ++m, ++m) {
				//(x,y)是16個種子點的位置
				float x = featcenter[m];
				float y = featcenter[m + 1];
				x_wght[m / 2] = max(1 - (fabs(x - x_sample)*1.0 / GridSpacing), 0);
				y_wght[m / 2] = max(1 - (fabs(y - y_sample)*1.0 / GridSpacing), 0);

			}
			for (int m = 0; m < 16; ++m)
				for (int n = 0; n < 8; ++n)
					pos_wght[m * 8 + n] = x_wght[m] * y_wght[m];
			free(x_wght);
			free(y_wght);

			//计算方向的加权，首先旋转梯度场到主方向，然后计算差异
			float diff[8], orient_wght[LEN];
			for (int m = 0; m < 8; ++m) {
				float angle = grad_sample - (p.ori) - orient_angles[m] + PI/**/; // 差值加上pi
				float temp = angle / (2.0 * PI/**/);
				angle -= (int)(temp)* (2.0 * PI/**/);
				diff[m] = angle - PI/**/;
			}

			// 计算高斯权重
			// Compute the gaussian weighting.
			float x = p.sx;
			float y = p.sy;
			float g = exp(-((x_sample - x)*(x_sample - x) + (y_sample - y)*(y_sample - y)) / (2 * feat_window*feat_window)) / (2 * PI/**/*feat_window*feat_window);

			// David G.Lowed的实验结果表明：对每个关键点，采用4*4*8共128维向量的描述子进项关键点表征，综合效果最佳
			// 计算对128维度的贡献值
			for (int m = 0; m < LEN; ++m) {
				//orient_wt是幅值方向在8個選定方向上的影響值。例如PI/3只對PI/4和PI/2有影響，權值視(0,PI/4)遠近從（1，0）漸變
				orient_wght[m] = max((1.0 - 1.0*fabs(diff[m % 8]) / orient_bin_spacing), 0);
				feat_desc[m] = feat_desc[m] + orient_wght[m] * pos_wght[m] * g*mag_sample;
			}
			free(pos_wght);
		}
		cout<<"get features"<<endl;
		free(feat);
		free(featcenter);
		//歸一化、抑制、再歸一化
		float norm = GetVectorNorm(feat_desc, LEN);
		for (int m = 0; m < LEN; m++) {
			feat_desc[m] /= norm;
			if (feat_desc[m] > 0.2)
				feat_desc[m] = 0.2;
		}
		norm = GetVectorNorm(feat_desc, LEN);
		for (int m = 0; m < LEN; m++) {
			feat_desc[m] /= norm;
			//printf("%f  ", feat_desc[m]);
		}
		//printf("\n");
		p.descrip = feat_desc;
	}
	free(feat_grid);
	free(feat_samples);
	cout<<"return"<<endl;
}

//为了显示图象金字塔，而作的图像水平拼接  
CImg<float> MySift::MosaicHorizen(CImg<float> im1, CImg<float> im2)
{
	int row, col;
	CImg<float> mosaic((im1.width() + im2.width()), max(im1.height(), im2.height()), 1, 1, 0);
	/* Copy images into mosaic1. */
	for (row = 0; row < im1.height(); row++)
		for (col = 0; col < im1.width(); col++)
			mosaic(col, row) = im1(col, row);
	for (row = 0; row < im2.height(); row++)
		for (col = 0; col < im2.width(); col++)
			mosaic((col + im1.width()), row) = im2(col, row);
	return mosaic;
}

//为了显示图象金字塔，而作的图像垂直拼接  
CImg<float> MySift::MosaicVertical(CImg<float> im1, CImg<float> im2) {
	int row, col;
	CImg<float> mosaic(max(im1.width(), im2.width()), im1.height() + im2.height(), 1, 1, 0);

	/* Copy images into mosaic1. */
	for (row = 0; row < im1.height(); row++)
		for (col = 0; col < im1.width(); col++)
			mosaic(col, row) = im1(col, row);
	for (row = 0; row < im2.height(); row++)
		for (col = 0; col < im2.width(); col++)
			mosaic(col, (row + im1.height())) = im2(col, row);
	return mosaic;
}

CImg<float> MySift::toGrayScale(CImg<float> img) {
    CImg<float> grayscaled(img.width(), img.height(), 1, 1, 0);
    cimg_forXY(img, x, y) {
        int r = img(x, y, 0);
        int g = img(x, y, 1);
        int b = img(x, y, 2);
        float newValue = (r * 0.2126 + g * 0.7152 + b * 0.0722);
        grayscaled(x, y) = newValue;
    }
    return grayscaled;
}

void MySift::cvConvertScale(CImg<float> src, CImg<float> &dst, float scale, float shift){
    cimg_forXY(src, x, y){
        dst(x, y) = src(x, y) * scale + shift;
    }
}
void MySift::cvConvertScaleAbs(CImg<float> src, CImg<float> &dst, float scale, float shift){
    cimg_forXY(src, x, y){
        dst(x, y) = abs(src(x, y) * scale + shift);
    }
}

void MySift::cvMinMaxLoc(CImg<float> img, float *min, float *max){
    float minTmp = FLT_MAX, maxTmp = FLT_MIN;
    cimg_forXY(img, x, y)
	{
		if (img(x, y) > maxTmp) maxTmp = img(x, y);
		if (img(x, y) < minTmp) minTmp = img(x, y);
	}
	*min = minTmp;
    *max = maxTmp;
}

void MySift::cvAddS(CImg<float> img, float input, CImg<float> &dst){
    cimg_forXY(dst, x, y){
        dst(x, y) = img(x, y) + input;
    }
}

void MySift::SiftMainProcess() {
	//读取图片
	CImg<float> tmp(filename);
	src = CImg<float>(tmp.width(), tmp.height(), 1, 3, 0);
	cimg_forXY(src, x, y)
	{
		src(x, y, 0) = (float) tmp(x, y, 0);
		src(x, y, 1) = (float) tmp(x, y, 1);
		src(x, y, 2) = (float) tmp(x, y, 2);
	}

	// //为图像分配内存   
	// image_kp = cvCreateImage(cvSize(src.width, src.height), IPL_DEPTH_8U, 3);
	// image_featDir = cvCreateImage(cvSize(src.width, src.height), IPL_DEPTH_8U, 3);
	// grey_im1 = cvCreateImage(cvSize(src.width, src.height), IPL_DEPTH_8U, 1);
	//doubleSizeImage = cvCreateImage(cvSize(2 * (src.width), 2 * (src.height)), IPL_DEPTH_8U, 3);

	//为图像阵列分配内存，假设两幅图像的大小相同，tempMat跟随image1的大小  
	//image1Mat = cvCreateMat(src.height, src.width, CV_32FC1);


	//转化成单通道图像再处理 从rgb转换到gray
	grey_im1 = toGrayScale(src);
	//转换进入Mat数据结构,图像操作使用的是浮点型操作  将图像数据转换为Mat格式的数据
	//cvConvert(grey_im1, image1Mat);
	CImg<float> image1Mat = grey_im1;
	//float t = (float)cvGetTickCount();


	//图像归一化  
	cvConvertScale(image1Mat, image1Mat, 1.0 / 255, 0);
	int dim = min(image1Mat.height(), image1Mat.width());
	numoctaves = (int) (log((float) dim) / log(2.0)) - 2;    //金字塔阶数
	numoctaves = min(numoctaves, MAXOCTAVES);

	//SIFT算法第一步，预滤波除噪声，建立金字塔底层  
	tempMat = ScaleInitImage(image1Mat);
	//SIFT算法第二步，建立Guassian金字塔和DOG金字塔
	BuildGaussianOctaves(tempMat);

	//t = (float)cvGetTickCount() - t;
	//printf("the time of build Gaussian pyramid and DOG pyramid is %.1f\n", t / (cvGetTickFrequency()*1000.));

	//显示高斯金字塔  
	for (int i = 0; i < numoctaves; i++) {
		if (i == 0) {
			mosaicHorizen1 = MosaicHorizen((Gaussianpyr[0].Octave)[0].Level, (Gaussianpyr[0].Octave)[1].Level);
			for (int j = 2; j < SCALESPEROCTAVE + 3; j++)
				mosaicHorizen1 = MosaicHorizen(mosaicHorizen1, (Gaussianpyr[0].Octave)[j].Level);
			for (int j = 0; j < NUMSIZE; j++)
				mosaicHorizen1 = halfSizeImage(mosaicHorizen1);
		} else if (i == 1) {
			mosaicHorizen2 = MosaicHorizen((Gaussianpyr[1].Octave)[0].Level, (Gaussianpyr[1].Octave)[1].Level);
			for (int j = 2; j < SCALESPEROCTAVE + 3; j++)
				mosaicHorizen2 = MosaicHorizen(mosaicHorizen2, (Gaussianpyr[1].Octave)[j].Level);
			for (int j = 0; j < NUMSIZE; j++)
				mosaicHorizen2 = halfSizeImage(mosaicHorizen2);
			mosaicVertical1 = MosaicVertical(mosaicHorizen1, mosaicHorizen2);
		} else {
			mosaicHorizen1 = MosaicHorizen((Gaussianpyr[i].Octave)[0].Level, (Gaussianpyr[i].Octave)[1].Level);
			for (int j = 2; j < SCALESPEROCTAVE + 3; j++)
				mosaicHorizen1 = MosaicHorizen(mosaicHorizen1, (Gaussianpyr[i].Octave)[j].Level);
			for (int j = 0; j < NUMSIZE; j++)
				mosaicHorizen1 = halfSizeImage(mosaicHorizen1);
			mosaicVertical1 = MosaicVertical(mosaicVertical1, mosaicHorizen1);
		}
	}
	mosaic1 = CImg<float>(mosaicVertical1.width(), mosaicVertical1.height(), 1, 1, 0);
	cvConvertScale(mosaicVertical1, mosaicVertical1, 255.0, 0);
	cvConvertScaleAbs(mosaicVertical1, mosaic1, 1, 0);
	mosaic1.display();


	//显示DOG金字塔  
	for (int i = 0; i < numoctaves; i++) {
		if (i == 0) {
			mosaicHorizen1 = MosaicHorizen((DOGoctaves[0].Octave)[0].Level, (DOGoctaves[0].Octave)[1].Level);
			for (int j = 2; j < SCALESPEROCTAVE + 2; j++)
				mosaicHorizen1 = MosaicHorizen(mosaicHorizen1, (DOGoctaves[0].Octave)[j].Level);
			for (int j = 0; j < NUMSIZE; j++)
				mosaicHorizen1 = halfSizeImage(mosaicHorizen1);
		} else if (i == 1) {
			mosaicHorizen2 = MosaicHorizen((DOGoctaves[1].Octave)[0].Level, (DOGoctaves[1].Octave)[1].Level);
			for (int j = 2; j < SCALESPEROCTAVE + 2; j++)
				mosaicHorizen2 = MosaicHorizen(mosaicHorizen2, (DOGoctaves[1].Octave)[j].Level);
			for (int j = 0; j < NUMSIZE; j++)
				mosaicHorizen2 = halfSizeImage(mosaicHorizen2);
			mosaicVertical1 = MosaicVertical(mosaicHorizen1, mosaicHorizen2);
		} else {
			mosaicHorizen1 = MosaicHorizen((DOGoctaves[i].Octave)[0].Level, (DOGoctaves[i].Octave)[1].Level);
			for (int j = 2; j < SCALESPEROCTAVE + 2; j++)
				mosaicHorizen1 = MosaicHorizen(mosaicHorizen1, (DOGoctaves[i].Octave)[j].Level);
			for (int j = 0; j < NUMSIZE; j++)
				mosaicHorizen1 = halfSizeImage(mosaicHorizen1);
			mosaicVertical1 = MosaicVertical(mosaicVertical1, mosaicHorizen1);
		}
	}
	//考虑到DOG金字塔各层图像都会有正负，所以，必须寻找最负的，以将所有图像抬高一个台阶去显示  
	float min_val = 0;
	float max_val = 0;
	cvMinMaxLoc(mosaicVertical1, &min_val, &max_val);
	//printf("%f\n", min_val);
	if (min_val < 0.0)
		//！
		cvAddS(mosaicVertical1, (-1.0) * min_val, mosaicVertical1);
	mosaic2 = CImg<float>(mosaicVertical1.width(), mosaicVertical1.height(), 1, 1, 0);
	cvConvertScale(mosaicVertical1, mosaicVertical1, 255.0 / (max_val - min_val), 0);
	cvConvertScaleAbs(mosaicVertical1, mosaic2, 1, 0);
	mosaic2.display();

	//SIFT算法第三步：特征点位置检测，最后确定特征点的位置  
	DetectKeypoint(numoctaves, Gaussianpyr);
	printf("the keypoints number are %d ;\n", keypoints.size());
	image_kp = src;
	DisplayKeypointLocation(image_kp, Gaussianpyr);
	image_kp.display();

	//SIFT算法第四步：计算高斯图像的梯度方向和幅值，计算各个特征点的主方向  
	ComputeGrad_DirecandMag(numoctaves, Gaussianpyr);
	AssignTheMainOrientation(numoctaves, Gaussianpyr, mag_pyr, grad_pyr);
	image_featDir = src;
	// cvCopy(src, image_featDir, NULL);
	DisplayOrientation(image_featDir, Gaussianpyr);

	//SIFT算法第五步：抽取各个特征点处的特征描述字  
	ExtractFeatureDescriptors(numoctaves, Gaussianpyr);
}

vector<Keypoint> MySift::getFirstKeyDescriptors() {
	return keyDescriptors;
}

void MySift::saveImgWithKeypoint(char* filename){
    image_kp.save(filename);
}