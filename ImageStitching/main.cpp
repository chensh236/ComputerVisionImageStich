#include "MyMatching.h"
#include "MyBlending.h"
#include "Projection.h"
#include <map>
#define MAX_IMG 2
#define RESIZE_SIZE 500.0
extern "C" {
#include "vl/generic.h"
#include "vl/sift.h"
};
typedef struct Image{
	CImg<float> resizeImg;
	CImg<unsigned char> src;
	CImg<unsigned char> grayDst;
    map<vector<float>, VlSiftKeypoint> features;
	vector<Keypoint> keypoints2;
} Image;

CImg<unsigned char> get_gray_image(const CImg<unsigned char> &srcImg) {
    if (srcImg.spectrum() == 1) {
        return srcImg;
    }
    int width = srcImg.width();
    int height = srcImg.height();
    int depth = srcImg.depth();
    CImg<unsigned char> grayImg(width, height, depth, 1);
    float r, g, b, gr;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            r = srcImg(i, j, 0);
            g = srcImg(i, j, 1);
            b = srcImg(i, j, 2);
            gr = 0.299*(r)+0.587*(g)+0.114*(b);
            grayImg(i, j) = gr;
        }
    }
    return grayImg;
}

vector<Image> imgs;
void processImg(const CImg<unsigned char>& src, int index, const char* filename){
	// 图片大小变换
	CImg<float> dst(src);
	float resize_factor;

	if (src.width() < src.height()) {
		resize_factor = RESIZE_SIZE / src.width();
	}
	else {
		resize_factor = RESIZE_SIZE / src.height();
	}

	if (resize_factor >= 1) {
		resize_factor = 1;
	}
	else {
		dst.resize(dst.width() * resize_factor, dst.height() * resize_factor, dst.depth(), dst.spectrum(), 3);
	}
	//dst.display("resize");
	// 覆盖的sift算法
	imgs[index].resizeImg = dst;
    imgs[index].grayDst = get_gray_image(imgs[index].resizeImg);
    imgs[index].grayDst.display();
	MySift mySift(filename, 1, dst);
	mySift.SiftMainProcess();
	imgs[index].keypoints2 = mySift.getFirstKeyDescriptors();
	for(Keypoint k : imgs[index].keypoints2){
		k.col /= resize_factor;
		k.row /= resize_factor;
		k.ix = k.col;
		k.iy = k.row;
	}

//    // 使用vl的sift算法
//    vl_sift_pix *imageData = new vl_sift_pix[dst.width()*dst.height()];
//
//    for (int i = 0; i < dst.width(); i++) {
//        for (int j = 0; j < dst.height(); j++) {
//            imageData[j*dst.width() + i] = dst(i, j, 0);
//        }
//    }
//
//    // Initialize a SIFT filter object.
//    int noctaves = 4, nlevels = 2, o_min = 0;
//    VlSiftFilt *siftFilt = NULL;
//    siftFilt = vl_sift_new(dst.width(), dst.height(), noctaves, nlevels, o_min);
//
//    map<vector<float>, VlSiftKeypoint> features;
//
//    // Compute the first octave of the DOG scale space.
//    if (vl_sift_process_first_octave(siftFilt, imageData) != VL_ERR_EOF) {
//        while (true) {
//            // Run the SIFT detector to get the keypoints.
//            vl_sift_detect(siftFilt);
//
//            VlSiftKeypoint *pKeyPoint = siftFilt->keys;
//            // For each keypoint:
//            for (int i = 0; i < siftFilt->nkeys; i++) {
//                VlSiftKeypoint tempKeyPoint = *pKeyPoint;
//
//                // Get the keypoint orientation(s).
//                double angles[4];
//                int angleCount = vl_sift_calc_keypoint_orientations(siftFilt, angles, &tempKeyPoint);
//
//                // For each orientation:
//                for (int j = 0; j < angleCount; j++) {
//                    double tempAngle = angles[j];
//                    vl_sift_pix descriptors[128];
//
//                    // Get the keypoint descriptor
//                    vl_sift_calc_keypoint_descriptor(siftFilt, descriptors, &tempKeyPoint, tempAngle);
//
//                    vector<float> des;
//                    int k = 0;
//                    while (k < 128) {
//                        des.push_back(descriptors[k]);
//                        k++;
//                    }
//
//                    tempKeyPoint.x /= resize_factor;
//                    tempKeyPoint.y /= resize_factor;
//                    tempKeyPoint.ix = tempKeyPoint.x;
//                    tempKeyPoint.iy = tempKeyPoint.y;
//
//                    features.insert(pair<vector<float>, VlSiftKeypoint>(des, tempKeyPoint));
//                }
//
//                pKeyPoint++;
//            }
//            if (vl_sift_process_next_octave(siftFilt) == VL_ERR_EOF) {
//                break;
//            }
//        }
//    }
//
//    vl_sift_delete(siftFilt);
//    delete[] imageData;
//    imageData = NULL;
//
//    for(auto iter = imgs[index].features.begin(); iter != imgs[index].features.end(); ++iter){
//        Keypoint newKey;
//        // keypoint转换
//        //newKey.descrip = (float*)malloc(128 * sizeof(float));
//        newKey.row = (iter->second).y;
//        newKey.col = (iter->second).y;
//        newKey.ix = (iter->second).ix;
//        newKey.iy = (iter->second).iy;
//        for(int i = 0; i < 128; i++){
//            newKey.descrip[i] = (iter->first)[i];
//        }
//        imgs[index].keypoints2.push_back(newKey);
//    }
}


int main() {
	for(int i = 0; i < MAX_IMG; i++){
		Image thisImg;
		imgs.push_back(thisImg);
		string inputFileName = "../ImageStitching/Input/" + to_string (i + 1) + ".bmp";
		imgs[i].src = CImg<unsigned char>(inputFileName.c_str());
		//求出投影
		CImg<unsigned char> projectedImg = Projection::imageProjection(imgs[i].src);
		processImg(projectedImg, i, inputFileName.c_str());
	}
    MyMatching myMatching(imgs[0].keypoints2,imgs[1].keypoints2);
    myMatching.featureMatchMainProcess();
    myMatching.mixImageAndDrawPairLine("../Output/mixImg.bmp", "../Output/mixImgWithLine.bmp");
    myMatching.drawOriKeypointOnImg(imgs[0].resizeImg, imgs[1].resizeImg, "../Output/1_kp_real.bmp", "../Output/2_kp_real.bmp");

//    for(int i = 0; i < MAX_IMG; ++i){
//        for(int j = 0; j < imgs[i].keypoints2.size(); ++j){
//            delete(imgs[i].keypoints2[j].descrip);
//            imgs[i].keypoints2[j].descrip = 0;
//        }
//    }

//	mySift2.saveImgWithKeypoint("../Output/2_kp.bmp");

	// MyMatching myMatching(mySift1.getFirstKeyDescriptors(),  mySift2.getFirstKeyDescriptors());
	// myMatching.featureMatchMainProcess();
	// myMatching.drawOriKeypointOnImg(inputAddr1, inputAddr2, "../Output/1_kp_real.bmp", "../Output/2_kp_real.bmp");
	// myMatching.mixImageAndDrawPairLine("../Output/mixImg.bmp", "../Output/mixImgWithLine.bmp");
	// myMatching.myRANSACtoFindKpTransAndDrawOut("../Output/mixImgWithLine_fixed.bmp");

	// MyBlending myBlending(myMatching.getMatchVec().col, myMatching.getMatchVec().row);
	// myBlending.blendingMainProcess(inputAddr1, inputAddr2);
	// myBlending.saveBlendedImg("../Output/blendedImg.bmp");

	return 0;
}
