#include <stdio.h>
#include <stdlib.h>
#include "gdal.h"
#include "cpl_conv.h"
#include "geotiffutil.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <regex>

using namespace std;

Geotiff readGeotiff(const char* fname) {

    GDALAllRegister();
    GDALDatasetH geotiffReader;
    GDALDriverH geotiffDriver;
    Geotiff GeotiffObj;

    GeotiffObj.filename = fname;
    geotiffReader = GDALOpen(fname, GA_ReadOnly);
    geotiffDriver = GDALGetDatasetDriver(geotiffReader);


    readDimensionsGeotiff(fname, geotiffReader, &GeotiffObj);
    readProjectionGeotiff(geotiffReader, &GeotiffObj);
    readBandGeotiff(fname, geotiffReader, &GeotiffObj, 1);


    return GeotiffObj;

}

void readBandGeotiff(
    const char* filename,
    GDALDatasetH reader,
    Geotiff* tiff,
    int bandNumber
) {

    int ncols, nrows;
    int row;
    int col;

    GDALRasterBandH handleBand;
    float** outBand;
    handleBand = GDALGetRasterBand(reader, bandNumber);

   
    ncols = tiff->xsize;
    nrows = tiff->ysize;

    outBand = (float**)malloc(nrows * sizeof(float*));
    GDALDataType datatype;
    datatype = GDALGetRasterDataType(GDALGetRasterBand(reader, bandNumber));
    cout << datatype << endl;
    float* scanLine = (float*) CPLMalloc(ncols* sizeof(float));

    for (row = 0; row < nrows; row++) {
        outBand[row] = (float*)malloc(ncols * sizeof(float));
        GDALRasterIO(
            handleBand, GF_Read, 0, row, ncols, 1,
            scanLine, ncols, 1,
            GDT_Float32, 0, 0
        );

        for (col = 0; col < ncols; col++) {
            outBand[row][col] = (float)scanLine[col];
        }
    }

    CPLFree(scanLine);
    tiff->band = outBand;

}


void readProjectionGeotiff(
    GDALDatasetH reader,
    Geotiff* tiff
) {

    double gt[DIMGT];
    double* pGT; 
    char* projection;
    GDALGetGeoTransform(reader, gt);
    pGT = &gt[0]; 
    tiff->projection = (char*)GDALGetProjectionRef(reader);

    int counter; 
    counter = 0;
    while (counter < DIMGT) {
        tiff->geotransform[counter] = *(pGT + counter);
        counter++;
    }
}



void readDimensionsGeotiff(
    const char* filename,
    GDALDatasetH reader,
    Geotiff* tiff
) {

    int xdim, ydim;
    int nbands;
    xdim = GDALGetRasterXSize(reader);
    ydim = GDALGetRasterYSize(reader);
    nbands = GDALGetRasterCount(reader);

    tiff->xsize = xdim;
    tiff->ysize = ydim;
    tiff->nbands = nbands;

}


GeotiffMeta readGeoTiffMeta(string& filename, int band) {
    ifstream file(filename);

    string band_num = to_string(band);

    string line;
    GeotiffMeta metaobj;

    string RADI_MUL_key = "RADIANCE_MULT_BAND_" + band_num;
    string RADI_ADD_key = "RADIANCE_ADD_BAND_" + band_num;
    string REFL_MUL_key = "REFLECTANCE_MULT_BAND_" + band_num;
    string REFL_ADD_key = "REFLECTANCE_ADD_BAND_" + band_num;
    string SUN_ELEV_key = "SUN_ELEVATION";
    string K1_CONSTANT_key = "K1_CONSTANT_BAND_" + band_num;
    string K2_CONSTANT_key = "K2_CONSTANT_BAND_" + band_num;

    float RADI_MUL =NULL;
    float RADI_ADD = NULL;
    float REFL_MUL = NULL;
    float REFL_ADD = NULL;
    float SUN_ELEV = NULL;
    float K1_CONSTANT = NULL;
    float K2_CONSTANT = NULL;

    if (band < 10) {
        while (getline(file, line)) {
            // 따옴표 제거
            line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
            // 각 키에 맞는 값 찾기
            if (line.find(RADI_MUL_key) != std::string::npos) {
                cout << "----------------------------------" << endl;
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    RADI_MUL = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(RADI_ADD_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    RADI_ADD = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(REFL_MUL_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    REFL_MUL = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(REFL_ADD_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    REFL_ADD = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(SUN_ELEV_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    SUN_ELEV = std::stof(line.substr(pos + 1));
                }
            }
        }
    }
    else if (band == 10 || band == 11) {
        while (getline(file, line)) {
            // 따옴표 제거
            line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
            // 각 키에 맞는 값 찾기
            if (line.find(RADI_MUL_key) != std::string::npos) {
                cout << "----------------------------------" << endl;
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    RADI_MUL = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(RADI_ADD_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    RADI_ADD = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(K1_CONSTANT_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    K1_CONSTANT = std::stof(line.substr(pos + 1));
                }
            }
            else if (line.find(K2_CONSTANT_key) != std::string::npos) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    K2_CONSTANT = std::stof(line.substr(pos + 1));
                }
            }
        }
    }


    const char* fname = filename.c_str();
    if (fname != NULL){
        metaobj.filename = fname;
    }
    else {
        metaobj.filename = NULL;
    }
    if (RADI_MUL != NULL) {
        metaobj.radi_mul = RADI_MUL;
    }
    else {
        metaobj.radi_mul = NULL;
    }
    if (RADI_ADD != NULL) {
        metaobj.radi_add = RADI_ADD;
    }
    else {
        metaobj.radi_add = NULL;
    }
    if (REFL_MUL != NULL) {
        metaobj.refl_mul = REFL_MUL;
    }
    else {
        metaobj.refl_mul = NULL;
    }
    if (REFL_ADD != NULL) {
        metaobj.refl_add = REFL_ADD;
    }
    else {
        metaobj.refl_add = NULL;
    }
    if (SUN_ELEV != NULL) {
        metaobj.sun_elev = SUN_ELEV;
    }
    else {
        metaobj.sun_elev = NULL;
    }
    if (K1_CONSTANT != NULL) {
        metaobj.k1_constant = K1_CONSTANT;
    }
    else {
        metaobj.k1_constant = NULL;
    }
    if (K2_CONSTANT != NULL) {
        metaobj.k2_constant = K2_CONSTANT;
    }
    else {
        metaobj.k2_constant = NULL;
    }
    return metaobj;
}


int writeGeotiff(
    Geotiff bandRef,
    vector<float**>bands,
    string& fname
) {

    char* projection = bandRef.projection;
    double geotransform[DIMGT];
    GDALDriverH outHandleDriver;
    GDALDatasetH outDataset;
    int nrows, ncols;
    int row, col;
    ncols = bandRef.xsize;
    nrows = bandRef.ysize;
    const char* outname = fname.c_str();

    int band_cnt = bands.size();

    int geo_count = 0;
    while (geo_count < DIMGT) {
        geotransform[geo_count] = bandRef.geotransform[geo_count];
        geo_count++;
    }

    outHandleDriver = GDALGetDriverByName("GTiff");
    outDataset = GDALCreate(outHandleDriver,
        outname,
        ncols, nrows, band_cnt,
        GDT_Float32, NULL);
    GDALSetGeoTransform(outDataset, geotransform);
    GDALSetProjection(outDataset, projection);


    float* scanLinebands[11];

    for (int i = 0; i < band_cnt; i++) {
        scanLinebands[i] = (float*)CPLMalloc(sizeof(float) * ncols);

    }

    GDALRasterBandH handleBands[11];


    for (int i = 0; i < band_cnt; i++) {
        handleBands[i] = GDALGetRasterBand(outDataset, i + 1);
    }
    for (int i = 0; i < band_cnt; i++) {
        handleBands[i] = GDALGetRasterBand(outDataset, i + 10);
    }

    for (row = 0; row < nrows; row++) {
        // 각 밴드의 데이터를 하나의 배열에 복사
        for (int i = 0; i < band_cnt; i++) {
            for (int col = 0; col < ncols; col++) {
                scanLinebands[i][col] = bands[i][row][col]; // 각 밴드에서 데이터 복사
            }
        }

        // 각 밴드를 한번에 처리
        for (int i = 0; i < band_cnt; i++) {
            GDALRasterIO(
                handleBands[i], GF_Write, 0, row, ncols, 1,
                scanLinebands[i], ncols, 1,
                GDT_Float32, 0, 0
            );
        }
    }


    GDALClose(outDataset);
    return 0;
}

int writeGeoband(
    Geotiff bandRef,
    float** band,
    string& fname
    ) {

    char* projection = bandRef.projection;
    double geotransform[DIMGT];
    GDALDriverH outHandleDriver;
    GDALDatasetH outDataset;
    int nrows, ncols;
    int row, col;
    ncols = bandRef.xsize;
    nrows = bandRef.ysize;
    const char* outname = fname.c_str();

    int band_cnt = 1;

    int geo_count = 0;
    while (geo_count < DIMGT) {
        geotransform[geo_count] = bandRef.geotransform[geo_count];
        geo_count++;
    }

    outHandleDriver = GDALGetDriverByName("GTiff");
    outDataset = GDALCreate(outHandleDriver,
        outname,
        ncols, nrows, band_cnt,
        GDT_Float32, NULL);
    GDALSetGeoTransform(outDataset, geotransform);
    GDALSetProjection(outDataset, projection);


    float* scanLineband;
    scanLineband = (float*)CPLMalloc(sizeof(float) * ncols);



    GDALRasterBandH handleBands;
    handleBands = GDALGetRasterBand(outDataset, 1);


    for (row = 0; row < nrows; row++) {
        // 각 밴드의 데이터를 하나의 배열에 복사
        for (int col = 0; col < ncols; col++) {
            scanLineband[col] = band[row][col]; // 각 밴드에서 데이터 복사
        }
        GDALRasterIO(handleBands, GF_Write, 0, row, ncols, 1, scanLineband, ncols, 1, GDT_Float32, 0, 0);
    }

    GDALClose(outDataset);
    return 0;
}


int cvtToTOAreflectance(const char* img_band, string& metadata, string& result_path) {

    Geotiff band;
    band = readGeotiff(img_band);


    int ncol = band.xsize;
    int nrow = band.ysize;

    float** DN = band.band;

    string fname = string(img_band);

    regex pattern("B(\\d+)");
    smatch match;
    int band_num;

    if (regex_search(fname, match, pattern)) {
        band_num = stoi(match[1]);
    }
    else {
        cout << "Band number not found" << endl;
    }


    GeotiffMeta meta;
    meta = readGeoTiffMeta(metadata, band_num);


    float** TOA_refl = (float**)malloc(nrow * sizeof(float*));
    float sun_elev;
    sun_elev = meta.sun_elev;

    
    for (int row = 0; row < nrow; row++) {
        TOA_refl[row] = (float*)malloc(ncol * sizeof(float));

        for (int col = 0; col < ncol; col++) {
            TOA_refl[row][col] = ((DN[row][col] * meta.refl_mul) + meta.refl_add) / sin(sun_elev);

        }
    }


    writeGeoband(band, TOA_refl, result_path);
    return 0;

}

int cvtToRadiance(const char* img_band, string& metadata, string& result_path) {

    Geotiff band;
    band = readGeotiff(img_band);
    int ncol = band.xsize;
    int nrow = band.ysize;

    float** DN = band.band;


    string fname = string(img_band);

    regex pattern("B(\\d+)");
    smatch match;
    int band_num;

    if (regex_search(fname, match, pattern)) {
        band_num = stoi(match[1]);
    }
    else {
        cout << "Band number not found" << endl;
    }

    GeotiffMeta meta;
    meta = readGeoTiffMeta(metadata, band_num);


    float** radi = (float**)malloc(nrow * sizeof(float*));
    float sun_elev;
    sun_elev = meta.sun_elev;


    for (int row = 0; row < nrow; row++) {
        radi[row] = (float*)malloc(ncol * sizeof(float));
        for (int col = 0; col < ncol; col++) {
            radi[row][col] = ((DN[row][col] * meta.radi_mul) + meta.radi_add);
        }
    };

    writeGeoband(band, radi, result_path);
    return 0;

}

int cvtToTOABT(const char* img_bands[], string& metadata, string& result_path) {

    Geotiff b10_band, b11_band;

    b10_band = readGeotiff(img_bands[0]);
    b11_band = readGeotiff(img_bands[1]);
    
    int ncol = b11_band.xsize;
    int nrow = b11_band.ysize;

    float** b10_DN = b10_band.band;
    float** b11_DN = b11_band.band;
    

    GeotiffMeta b10_meta, b11_meta;
    for (int i = 0; i < 4; i++) {
        if (i == 0) {

            b10_meta = readGeoTiffMeta(metadata, 10);

        }
        else if (i == 1) {

            b11_meta = readGeoTiffMeta(metadata, 11);

        }
    }


    float** b10_thermal = (float**)malloc(nrow * sizeof(float*));
    float** b11_thermal = (float**)malloc(nrow * sizeof(float*));


    for (int row = 0; row < nrow; row++) {
        b10_thermal[row] = (float*)malloc(ncol * sizeof(float));
        b11_thermal[row] = (float*)malloc(ncol * sizeof(float));
       
        for (int col = 0; col < ncol; col++) {
            b10_thermal[row][col] = b10_meta.k2_constant/(log((b10_meta.k1_constant/((b10_DN[row][col] * b10_meta.radi_mul) + b10_meta.radi_add))+1));
            b11_thermal[row][col] = b11_meta.k2_constant / (log((b11_meta.k1_constant / ((b11_DN[row][col] * b11_meta.radi_mul) + b11_meta.radi_add)) + 1));

        }
    }

    vector<float**> outputbands;
    outputbands.push_back(b10_thermal);
    outputbands.push_back(b11_thermal);

    vector<string> tasks;
    tasks.push_back("radiance");
    tasks.push_back("reflectance");
    tasks.push_back("BT");

    writeGeotiff(b10_band,outputbands,result_path);
    return 0;

}