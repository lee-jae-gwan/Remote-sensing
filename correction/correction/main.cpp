#include <iostream>
#include "gdal_priv.h"
#include "GeotiffUtil.h"
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

int main() {
	GDALAllRegister();
	
	const char* img_bands[4];
	img_bands[0] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B2.TIF";
	img_bands[1] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B3.TIF";
	img_bands[2] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B4.TIF";
	img_bands[3] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B5.TIF";
	string metadata_path = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1\\/LC09_L1TP_115036_20250222_20250222_02_T1_MTL.txt";

	const char* thermal_bands[2];
	thermal_bands[0] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B10.TIF";
	thermal_bands[1] = "D:\\data\\LC09_L1TP_115036_20250222_20250222_02_T1/LC09_L1TP_115036_20250222_20250222_02_T1_B11.TIF";


	string result_path = "D:\\data\\correction\\result_img\\TOA_reflectance_bandsetTest2.TIF";
	string result_path2 = "D:\\data\\correction\\result_img\\Radiance_bandsetTest2.TIF";
	string result_path3 = "D:\\data\\correction\\result_img\\ThermalTemperature_bandsetTest2.TIF";




	cvtToRadiance(img_bands, metadata_path, result_path2);
	cvtToTOAreflectance(img_bands, metadata_path, result_path);
	cvtToTOABT(thermal_bands, metadata_path, result_path3);

}

