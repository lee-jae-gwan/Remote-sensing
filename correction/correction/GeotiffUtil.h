#pragma once
#ifndef GEOTIFFUTIL_H
#define GEOTIFFUTIL_H
#define DIMGT 6
#include<vector>

using namespace std; 

typedef struct {
	const char* filename;
	int xsize;
	int ysize;
	int nbands;
	double geotransform[6];
	double NoDataValue;
	char* projection;
	float** band;
}Geotiff;


typedef struct {
	const char* filename;
	float radi_mul;
	float radi_add;
	float refl_mul;
	float refl_add;
	float sun_elev;
	float k1_constant;
	float k2_constant;
} GeotiffMeta;


Geotiff readGeotiff(const char*);
void readProjectionGeotiff(GDALDatasetH, Geotiff*);
void readDimensionsGeotiff(const char*, GDALDatasetH, Geotiff*);
void readBandGeotiff(const char*, GDALDatasetH, Geotiff*, int);
GeotiffMeta readGeoTiffMeta(string&, int);
void geobandWrite(uint16_t**, Geotiff, string&);
int geotiffWrite(uint16_t**, uint16_t**, uint16_t**, uint16_t**, Geotiff, string&);

int writeGeotiff(
	Geotiff,
	vector<float**>,
	string&,
	string task);

int cvtToTOAreflectance(const char*[], string&, string&);
int cvtToRadiance(const char* [], string&, string&);
int cvtToTOABT(const char* [], string&, string&);

#endif