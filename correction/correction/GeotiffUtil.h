#pragma once
#ifndef GEOTIFFUTIL_H
#define GEOTIFFUTIL_H
#define DIMGT 6

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
	float**,
	float**,
	float**,
	float**,
	string&);

int TOA_reflectance(const char*[], string&, string&);

#endif