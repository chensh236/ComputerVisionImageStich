void DealWithImgData(BYTE *srcdata, BYTE *drcdata,int width,int height)//参数一为原图像的数据区首指针，参数二为投影后图像的数据区首指针，参数三为图像的宽，参数四为图像的高
{
	//双线性插值算法
	int i_original_img_hnum, i_original_img_wnum;//目标点坐标
	double distance_to_a_y, distance_to_a_x;//在原图像中与a点的水平距离  
	int original_point_a, original_point_b, original_point_c, original_point_d;
 
	int l_width = WIDTHBYTES(width* 24);//计算位图的实际宽度并确保它为4byte的倍数
	int drcpoint;
	double R = 1200;//像素距离
	double x, y;
	for (int hnum = 0; hnum < height; hnum++)
	{
		for (int wnum = 0; wnum < width; wnum++)
		{
			drcpoint = l_width*hnum + wnum * 3;//数组位置偏移量，对应于图像的各像素点RGB的起点
			//柱面投影
			double k = R / sqrt(R*R + (wnum- width / 2) * (wnum - width / 2));
			x = (wnum - width / 2) / k + width / 2;
			y = (hnum - height / 2) / k + height / 2;
			if (x >= 0 && y >= 0 && x < width && y < height)
			{
				/***********双线性插值算法***********/
				i_original_img_hnum = y;
				i_original_img_wnum = x;
				distance_to_a_y = y - i_original_img_hnum;
				distance_to_a_x = x - i_original_img_wnum;//在原图像中与a点的垂直距离  
 
				original_point_a = i_original_img_hnum*l_width + i_original_img_wnum * 3;//数组位置偏移量，对应于图像的各像素点RGB的起点,相当于点A    
				original_point_b = original_point_a + 3;//数组位置偏移量，对应于图像的各像素点RGB的起点,相当于点B  
				original_point_c = original_point_a + l_width;//数组位置偏移量，对应于图像的各像素点RGB的起点,相当于点C   
				original_point_d = original_point_c + 3;//数组位置偏移量，对应于图像的各像素点RGB的起点,相当于点D  
 
				if (i_original_img_hnum == height - 1)
				{
					original_point_c = original_point_a;
					original_point_d = original_point_b;
				}
				if (i_original_img_wnum == width - 1)
				{
					original_point_a = original_point_b;
					original_point_c = original_point_d;
				}
 
				drcdata[drcpoint + 0] =
					srcdata[original_point_a + 0] * (1 - distance_to_a_x)*(1 - distance_to_a_y) +
					srcdata[original_point_b + 0] * distance_to_a_x*(1 - distance_to_a_y) +
					srcdata[original_point_c + 0] * distance_to_a_y*(1 - distance_to_a_x) +
					srcdata[original_point_c + 0] * distance_to_a_y*distance_to_a_x;
				drcdata[drcpoint + 1] =
					srcdata[original_point_a + 1] * (1 - distance_to_a_x)*(1 - distance_to_a_y) +
					srcdata[original_point_b + 1] * distance_to_a_x*(1 - distance_to_a_y) +
					srcdata[original_point_c + 1] * distance_to_a_y*(1 - distance_to_a_x) +
					srcdata[original_point_c + 1] * distance_to_a_y*distance_to_a_x;
				drcdata[drcpoint + 2] =
					srcdata[original_point_a + 2] * (1 - distance_to_a_x)*(1 - distance_to_a_y) +
					srcdata[original_point_b + 2] * distance_to_a_x*(1 - distance_to_a_y) +
					srcdata[original_point_c + 2] * distance_to_a_y*(1 - distance_to_a_x) +
					srcdata[original_point_c + 2] * distance_to_a_y*distance_to_a_x;
				/***********双线性插值算法***********/
			}
		}
	}
}