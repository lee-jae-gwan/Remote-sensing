#include <iostream>
#include "gdal_priv.h"
#include "GeotiffUtil.h"
#include <fstream>
#include <string>
#include <cmath>
#include <regex>

using namespace std;

int main() {
	GDALAllRegister();

	const char* img_bands[9];
	img_bands[0] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B1.TIF";
	img_bands[1] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B2.TIF";
	img_bands[2] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B3.TIF";
	img_bands[3] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B4.TIF";
	img_bands[4] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B5.TIF";
	img_bands[5] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B6.TIF";
	img_bands[6] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B7.TIF";
	img_bands[7] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B8.TIF";
	img_bands[8] = "D:\\data\\LC09_L1TP_115036_20250121_20250121_02_T1/LC09_L1TP_115036_20250121_20250121_02_T1_B9.TIF";

	string metadata_path = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1\\/LC09_L1TP_115036_20250222_20250222_02_T1_MTL.txt";

	const char* thermal_bands[2];
	thermal_bands[0] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B10.TIF";
	thermal_bands[1] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B11.TIF";


	string result_path = "D:\\data\\correction\\result_img\\202501/";


	for (int i = 0; i < 9; i++) {
		string fname = string(img_bands[i]);
		string tmp = fname.substr(49);
		string tmp_path = result_path + "TOAreflectance_" + tmp;
		string tmp_path2 = result_path + "TOARadiance_" + tmp;

		cvtToTOAreflectance(img_bands[i], metadata_path, tmp_path);
		cvtToRadiance(img_bands[i], metadata_path, tmp_path2);
	}


	//cvtToTOABT(thermal_bands, metadata_path, result_path3);

}

