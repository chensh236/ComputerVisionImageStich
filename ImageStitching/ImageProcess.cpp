#include "ImageProcess.h"

// 拼合过程
ImageProcess::ImageProcess(string fileDic, const int picSum){
    // 文件处理与sift提取
    readFile(fileDic, picSum);
	matching();
    
}

// 文件读取
void ImageProcess::readFile(string fileDic, const int picSum){
    for(int i = 0; i < picSum; i++){
		Image currentImg;
        // format: fileDic/index.bmp
        // 顺序读取
		string inputFileName = fileDic + to_string (i + 1) + ".bmp";
		//求出投影
		currentImg.projectedSrc = Projection::imageProjection(CImg<unsigned char>(inputFileName.c_str()));
        // 转换为灰度图
        currentImg.graySrc = toGrayScale(currentImg.projectedSrc);
        currentImg.resizeImg = currentImg.graySrc;
        imgs.push_back(currentImg);
        // sift计算特征
        imgs[i].features = siftAlgorithm(i);
	}
}

// 转换灰度图像
CImg<unsigned char> ImageProcess::toGrayScale(const CImg<unsigned char>& src) {
    // 如果是灰度图像返回
    if (src.spectrum() == 1) {
        return src;
    }
    // 创建灰度图像
    CImg<unsigned char> gray(src.width(), src.height(), src.depth(), 1, 0);
    // 灰度图转换

    cimg_forXY(src, x, y){
        gray(x, y) = 0.299*(float)src(x, y, 0)+0.587*(float)src(x, y, 1)+0.114*(float)src(x, y, 2);
    }

    return gray;
}

// sift特征提取
//https://www.cnblogs.com/wangguchangqing/p/9176103.html
map<vector<float>, VlSiftKeypoint> ImageProcess::siftAlgorithm(const int index) {

    // 判断高度和宽度的比例
    bool flag = (imgs[index].graySrc.width() < imgs[index].graySrc.height());

    // 尺寸变化的比例
	float rf = flag ? (RESIZE_FACTOR / imgs[index].graySrc.width()) : 
			   (RESIZE_FACTOR / imgs[index].graySrc.height());

	if (rf >= 1){
		rf = 1;
	}
	else{
		imgs[index].resizeImg.resize(imgs[index].resizeImg.width() * rf, imgs[index].resizeImg.height() * rf, 
            imgs[index].resizeImg.depth(), imgs[index].resizeImg.spectrum(), 3);
	}

    // 新建vlFeat的图像
	vl_sift_pix *vlImg = new vl_sift_pix[imgs[index].resizeImg.width() * imgs[index].resizeImg.height()];

	for(int x = 0; x < imgs[index].resizeImg.width(); ++x){
		for(int y = 0; y < imgs[index].resizeImg.height(); ++y){
			vlImg[y * imgs[index].resizeImg.width() + x] = (float)imgs[index].resizeImg(x, y, 0);
		}
	}

	// Initialize a SIFT filter object.
	int noctaves = NOTAVES_NUM, nlevels = LEVEL_NUM, o_min = 0;
	VlSiftFilt *siftFilt = NULL;
	siftFilt = vl_sift_new(imgs[index].resizeImg.width(), imgs[index].resizeImg.height(), noctaves, nlevels, o_min);

	map<vector<float>, VlSiftKeypoint> features;

	// Compute the first octave of the DOG scale space.
	if (vl_sift_process_first_octave(siftFilt, vlImg) != VL_ERR_EOF) {
		while (true) {
			// Run the SIFT detector to get the keypoints.
			vl_sift_detect(siftFilt);

            // 通过高斯金字塔获取关键点
			VlSiftKeypoint *pKeyPoint = siftFilt->keys;
			// For each keypoint:
			for (int i = 0; i < siftFilt->nkeys; i++) {
				VlSiftKeypoint tempKeyPoint = *pKeyPoint;

				// Get the keypoint orientation(s).
				double angles[4];
				int angleCount = vl_sift_calc_keypoint_orientations(siftFilt, angles, &tempKeyPoint);

				// For each orientation:
				for (int j = 0; j < angleCount; j++) {
					double tempAngle = angles[j];
                    // 初始化描述子
					vl_sift_pix descriptors[128];

					// Get the keypoint descriptor
					vl_sift_calc_keypoint_descriptor(siftFilt, descriptors, &tempKeyPoint, tempAngle);

					vector<float> descriptorVec;
					int k = 0;
					while (k < 128) {
						descriptorVec.push_back(descriptors[k]);
						k++;
					}

					tempKeyPoint.x /= rf;
					tempKeyPoint.y /= rf;
					tempKeyPoint.ix = tempKeyPoint.x;
					tempKeyPoint.iy = tempKeyPoint.y;

					features.insert(pair<vector<float>, VlSiftKeypoint>(descriptorVec, tempKeyPoint));
				}

				pKeyPoint++;
			}
			if (vl_sift_process_next_octave(siftFilt) == VL_ERR_EOF) {
				break;
			}
		}
	}
    // 销毁指针
	vl_sift_delete(siftFilt);
	delete[] vlImg;
	vlImg = NULL;
	return features;
}

void ImageProcess::matching(){
	// 用于判断两张图像是否需要进行拼接
	stichingMat = new bool* [imgs.size()];
	for(int i = 0; i < imgs.size(); i++){
		stichingMat[i] = new bool[imgs.size()];
	}

	vector< vector<int> > featureMatchingMat(imgs.size());

	// 对所有图像进行遍历
	for(int i = 0; i < imgs.size(); ++i){
		for(int j = 0; j < imgs.size(); ++j){
			// 初始化
			stitchingMatp[i][j] = false;
			featureMatchingMat[i][j] = 0;

			if (i == j)
				continue;

			vector<pair> pairs = getPair(imgs[i], imgs[j]);
			if (pairs.size() >= THRESHOLD) {
				stichingMat[i][j] = true;
				matching_index[i].push_back(j);
			}
		}
	}
}

vector<pair> getPair(const Image& imgA, const Image& imgB){
	
	// 这里使用KD树节省查找时间
// dataType	type of data (VL_TYPE_FLOAT or VL_TYPE_DOUBLE)
// dimension	data dimensionality.
// numTrees	number of trees in the forest.
// distance	type of distance norm (VlDistanceL1 or VlDistanceL2). 这里选择L1范数，但是我在另外一篇看他使用的是L2范数

	VlKDForest* forest = vl_kdforest_new(VL_TYPE_FLOAT, 128, 1, VlDistanceL1);

	// 共有128 * number of keypoints 的特征点
	float *data = new float[128 * imgA.features.size()];

	int k = 0;
	// 将其放入到数组中
	for (auto it = imgA.features.begin(); it != imgA.features.end(); it++) {
		const vector<float> &descriptors = it->first;
		for (int i = 0; i < 128; i++) {
			data[128 * k + i] = descriptors[i];
		}
		k++;
	}

// Parameters
// self	KDTree object
// numData	number of data points.
// data	pointer to the data.
// The function builds the KDTree by processing the data data. 
// For efficiency, KDTree does not make a copy the data, but retains a pointer to it. 
// Therefore the data buffer must be valid and unchanged for the lifespan of the object.
// The number of data points numData must not be smaller than one.
	// 建立KD树
	vl_kdforest_build(forest, imgA.features.size(), data);

	vector<pair> res;

	// 利用KD树对特征点进行寻找
	VlKDForestSearcher* searcher = vl_kdforest_new_searcher(forest);
	// VlKDForestNeighbor : Neighbor of a query point.
	VlKDForestNeighbor neighbours[2];

	// 图像B的特征点数据
	for (auto it = feature_b.begin(); it != feature_b.end(); it++){
		float *temp_data = new float[128];

		for (int i = 0; i < 128; i++) {
			temp_data[i] = (it->first)[i];
		}

// Parameters
// self	object.
// neighbors	list of nearest neighbors found (output).
// numNeighbors	number of nearest neighbors to find.
// query	query point.
// Returns
// number of tree leaves visited.
// A neighbor is represented by an instance of the structure VlKDForestNeighbor. 
// Each entry contains the index of the neighbor (this is an index into the KDTree data) 
// and its distance to the query point. Neighbors are sorted by increasing distance.
		int nvisited = vl_kdforestsearcher_query(searcher, neighbours, 2, temp_data);

		float ratio = neighbours[0].distance / neighbours[1].distance;
		// 寻找到的两个点的距离比较远
		if (ratio < 0.5) {
			vector<float> des(128);
			// 利用A图像的128个描述子来反过来求出A的坐标
			for (int j = 0; j < 128; j++) {
				des[j] = data[j + neighbours[0].index * 128];
			}
			VlSiftKeypoint left = imgA.features.find(des)->second;
			VlSiftKeypoint right = it->second;
			// 获得匹配的对
			res.push_back(pair(left, right));
		}

		delete[] temp_data;
		temp_data = NULL;
	}

	vl_kdforestsearcher_delete(searcher);
	vl_kdforest_delete(forest);

	delete[] data;
	data = NULL;
	return res;
}