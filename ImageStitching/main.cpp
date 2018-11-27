#include "ImageProcess.h"

int main() {
	// for(int i = 0; i < MAX_IMG; i++){
	// 	Image thisImg;
	// 	imgs.push_back(thisImg);
	// 	string inputFileName = "../ImageStitching/Input/" + to_string (i + 1) + ".bmp";
	// 	imgs[i].src = CImg<unsigned char>(inputFileName.c_str());
	// 	//求出投影
	// 	CImg<unsigned char> projectedImg = Projection::imageProjection(imgs[i].src);
	// 	processImg(projectedImg, i, inputFileName.c_str());
	// }
    // MyMatching myMatching(imgs[0].keypoints2,imgs[1].keypoints2);
    // myMatching.featureMatchMainProcess();
    // myMatching.mixImageAndDrawPairLine("../Output/mixImg.bmp", "../Output/mixImgWithLine.bmp");
    // myMatching.drawOriKeypointOnImg(imgs[0].resizeImg, imgs[1].resizeImg, "../Output/1_kp_real.bmp", "../Output/2_kp_real.bmp");

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
	string ss = "../ImageStiching/Input/";
	ImageProcess ip(ss,2);
	return 0;
}
