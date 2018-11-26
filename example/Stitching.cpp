#include "stdafx.h"
#include "Stitching.h"
#include "Projection.h"
#include <queue>

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
			r = srcImg(i, j, 0, 0);
			g = srcImg(i, j, 0, 1);
			b = srcImg(i, j, 0, 2);
			gr = 0.299*(r)+0.587*(g)+0.114*(b);
			grayImg(i, j, 0, 0) = gr;
		}
	}
	return grayImg;
}

int getMiddleIndex(vector<vector<int>> matching_index) {
	int one_side = 0;

	for (int i = 0; i < matching_index.size(); i++) {
		if (matching_index[i].size() == 1) {
			one_side = i;
			break;
		}
	}

	int middle_index = one_side;
	int pre_middle_index = -1;
	int n = matching_index.size() / 2;
	
	while (n--) {
		for (int i = 0; i < matching_index[middle_index].size(); i++) {
			if (matching_index[middle_index][i] != pre_middle_index) {
				pre_middle_index = middle_index;
				middle_index = matching_index[middle_index][i];

				break;
			}
		}
	}

	return middle_index;
}

void updateFeaturesByHomography(map<vector<float>, VlSiftKeypoint> &feature, Parameters H, float offset_x, float offset_y) {
	for (auto iter = feature.begin(); iter != feature.end(); iter++) {
		float cur_x = iter->second.x;
		float cur_y = iter->second.y;
		iter->second.x = getXAfterWarping(cur_x, cur_y, H) - offset_x;
		iter->second.y = getYAfterWarping(cur_x, cur_y, H) - offset_y;
		iter->second.ix = int(iter->second.x);
		iter->second.iy = int(iter->second.y);
	}
}

void updateFeaturesByOffset(map<vector<float>, VlSiftKeypoint> &feature, int offset_x, int offset_y) {
	for (auto iter = feature.begin(); iter != feature.end(); iter++) {
		iter->second.x -= offset_x;
		iter->second.y -= offset_y;
		iter->second.ix = int(iter->second.x);
		iter->second.iy = int(iter->second.y);
	}
}

CImg<unsigned char> stitching(vector<CImg<unsigned char>> &src_imgs) {
	// Used to save each image's features and corresponding coordinates.
	vector<map<vector<float>, VlSiftKeypoint>> features(src_imgs.size());

	for (int i = 0; i < src_imgs.size(); i++) {
		cout << "Preprocessing input image " << i << " ..." << endl;

		cout << "Cylinder projection started." << endl;
		src_imgs[i] = cylinderProjection(src_imgs[i]);
		cout << "Cylinder projection finished." << endl;

		CImg<unsigned char> gray = get_gray_image(src_imgs[i]);

		cout << "Extracting SIFT feature started." << endl;
		features[i] = getFeatureFromImage(gray);
		cout << "Extracting SIFT feature finished." << endl;

		cout << "Preprocessing input image " << i << " finished." << endl << endl;
	}

	// Used to record the image's adjacent images.
	bool need_stitching[MAX_STITCHING_NUM][MAX_STITCHING_NUM] = { false };

	// Used to record each image's adjacent images.
	vector<vector<int>> matching_index(src_imgs.size());

	cout << "Finding adjacent images ..." << endl;

	for (int i = 0; i < src_imgs.size(); i++) {
		for (int j = 0; j < src_imgs.size(); j++) {
			if (i == j)
				continue;

			vector<point_pair> pairs = getPointPairsFromFeature(features[i], features[j]);
			if (pairs.size() >= 20) {
				need_stitching[i][j] = true;

				cout << "Image " << i << " and " << j << " are adjacent." << endl;
				matching_index[i].push_back(j);
			}

		}
	}

	cout << endl;

	cout << "Stitching begins." << endl << endl;
	// Stitching begins.

	// Stitching from middle to have better visual effect.
	int start_index = getMiddleIndex(matching_index);
	//int start_index = 3;

	// Used to record the previous stitched image.
	int prev_dst_index = start_index;

	queue<int> unstitched_index;
	unstitched_index.push(start_index);

	CImg<unsigned char> cur_stitched_img = src_imgs[start_index];

	while (!unstitched_index.empty()) {
		int src_index = unstitched_index.front();
		unstitched_index.pop();

		for (int j = matching_index[src_index].size() - 1; j >= 0; j--) {
			int dst_index = matching_index[src_index][j];

			if (need_stitching[src_index][dst_index] == false) {
				continue;
			}
			else {
				need_stitching[src_index][dst_index] = false;
				need_stitching[dst_index][src_index] = false;
				unstitched_index.push(dst_index);
			}

			// Matching features using best-bin kd-tree.
			vector<point_pair> src_to_dst_pairs = getPointPairsFromFeature(features[src_index], features[dst_index]);
			vector<point_pair> dst_to_src_pairs = getPointPairsFromFeature(features[dst_index], features[src_index]);

			if (src_to_dst_pairs.size() > dst_to_src_pairs.size()) {
				dst_to_src_pairs.clear();
				for (int i = 0; i < src_to_dst_pairs.size(); i++) {
					point_pair temp(src_to_dst_pairs[i].b, src_to_dst_pairs[i].a);
					dst_to_src_pairs.push_back(temp);
				}
			}
			else {
				src_to_dst_pairs.clear();
				for (int i = 0; i < dst_to_src_pairs.size(); i++) {
					point_pair temp(dst_to_src_pairs[i].b, dst_to_src_pairs[i].a);
					src_to_dst_pairs.push_back(temp);
				}
			}

			// Finding homography by RANSAC.
			Parameters forward_H = RANSAC(dst_to_src_pairs);
			Parameters backward_H = RANSAC(src_to_dst_pairs);

			// Calculate the size of the image after stitching.
			float min_x = getMinXAfterWarping(src_imgs[dst_index], forward_H);
			min_x = (min_x < 0) ? min_x : 0;
			float min_y = getMinYAfterWarping(src_imgs[dst_index], forward_H);
			min_y = (min_y < 0) ? min_y : 0;
			float max_x = getMaxXAfterWarping(src_imgs[dst_index], forward_H);
			max_x = (max_x >= cur_stitched_img.width()) ? max_x : cur_stitched_img.width();
			float max_y = getMaxYAferWarping(src_imgs[dst_index], forward_H);
			max_y = (max_y >= cur_stitched_img.height()) ? max_y : cur_stitched_img.height();

			int new_width = ceil(max_x - min_x);
			int new_height = ceil(max_y - min_y);

			CImg<unsigned char> a(new_width, new_height, 1, src_imgs[dst_index].spectrum(), 0);
			CImg<unsigned char> b(new_width, new_height, 1, src_imgs[dst_index].spectrum(), 0);

			// Warping the dst image into the coordinate space of currently stitched image.
			warpingImageByHomography(src_imgs[dst_index], a, backward_H, min_x, min_y);
			// Move stitched image according to min_x and min_y as translation.
			movingImageByOffset(cur_stitched_img, b, min_x, min_y);

			// Update features coordinates according to homography.
			updateFeaturesByHomography(features[dst_index], forward_H, min_x, min_y);
			// Update features coordinates according to min_x and min_y.
			updateFeaturesByOffset(features[prev_dst_index], min_x, min_y);

			// Blending two images.
			cur_stitched_img = blendTwoImages(a, b);
			prev_dst_index = dst_index;

			cout << "Image " << dst_index << " has stitched." << endl << endl;
		}
	}

	return cur_stitched_img;
}