#ifndef IMAGEPROCESS_H
#define IMAGEPROCESS_H
#include "CImg.h"
// 用于图像变形
#define RESIZE_FACTOR 500.0
// 用于sift
#define NOTAVES_NUM 4
#define LEVEL_NUM 2
// 用于图像匹配
#define THRESHOLD 20
#include <map>
#include <iostream>
#include <cstring>
#include <vector>
#include "Projection.h"

extern "C" {
    #include "vl/generic.h"
    #include "vl/sift.h"
};
using namespace std;
using namespace cimg_library;

// 用于存储关键点的对
struct pair {
	VlSiftKeypoint a;
	VlSiftKeypoint b;
	pair(VlSiftKeypoint _a, VlSiftKeypoint _b):a(_a), b(_b){}
};

// 存储图片的信息，比如cimg 关键点等
typedef struct Image{
    // 投影后的图像
    CImg<unsigned char> projectedSrc;
    // 灰度图像
	CImg<unsigned char> graySrc;
    // 尺寸转变后的图像
    CImg<unsigned char> resizeImg;
    // 存储关键点的特征描述子
    map<vector<float>, VlSiftKeypoint> features;
} Image;



// 这个类进行图像处理，包括读取图像、进行sift特征提取等
class ImageProcess{
public:
    ImageProcess(string, const int);
    ~ImageProcess(){
        // 销毁指针
        for(int i = 0; i < imgs.size(); i++){
            delete(stichingMat[i]);
            stichingMat[i] = 0;
        }
        delete(stichingMat);
        stichingMat = 0;
    }

    // 这个类进行文件读取操作
    // string : filename
    // int : sum of the pictures
    void readFile(string, const int);

private:
    vector<Image> imgs;
    CImg<unsigned char> toGrayScale(const CImg<unsigned char>&);
    // 利用sift获取特征点
    // int : 图片的序号
    map<vector<float>, VlSiftKeypoint> siftAlgorithm(const int);
    // 进行图像匹配
    void matching();
    // 求出匹配的关键点
    vector<pair> getPair(const Image&, const Image&);
    // 判断是否需要拼接
    bool **stichingMat
};

#endif