#include <stdio.h>
#include <stdlib.h>
#include "gdal.h"
#include "cpl_conv.h"
#include "geotiffutil.h"
#include <string>
#include <fstream>
#include <iostream>

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

    float RADI_MUL;
    float RADI_ADD;
    float REFL_MUL;
    float REFL_ADD;
    float SUN_ELEV;

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

    const char* fname = filename.c_str();
    metaobj.filename = fname;
    metaobj.radi_mul = RADI_MUL;
    metaobj.radi_add = RADI_ADD;
    metaobj.refl_mul = REFL_MUL;
    metaobj.refl_add = REFL_ADD;
    metaobj.sun_elev = SUN_ELEV;// 기본값
    
    return metaobj;
}

void geobandWrite(float** band, Geotiff src ,string& fname) {

    char* projection = src.projection;
    double geotransform[DIMGT];
    GDALDriverH outHandleDriver;
    GDALDatasetH outDataset;
    const char* outname = fname.c_str();

    int nrows, ncols;
    int row, col;
    ncols = src.xsize;
    nrows = src.ysize;

    int count;
    count = 0;
    while (count < DIMGT) {
        geotransform[count] = src.geotransform[count];
        count++;
    }
    outHandleDriver = GDALGetDriverByName("GTiff");
    outDataset = GDALCreate(outHandleDriver, outname, ncols, nrows, 1, GDT_UInt16, NULL);
    GDALSetGeoTransform(outDataset, geotransform);
    GDALSetProjection(outDataset, projection);

    uint16_t* scanLine = (uint16_t*)CPLMalloc(sizeof(uint16_t) * ncols);

    GDALRasterBandH handleRedBand;

    handleRedBand = GDALGetRasterBand(outDataset, 1);

    for (row = 0; row < nrows; row++) {
        for (col = 0; col < ncols; col++) {
            scanLine[col] = band[row][col];

        }
        GDALRasterIO(handleRedBand, GF_Write, 0, row, ncols, 1, scanLine, ncols, 1, GDT_UInt16, 0, 0);
    }
    GDALClose(outDataset);
}


int geotiffWrite(uint16_t** bband, uint16_t** gband, uint16_t** rband, uint16_t** nband, Geotiff src, string& fname) {

    char* projection = src.projection;
    double geotransform[DIMGT];
    GDALDriverH outHandleDriver;
    GDALDatasetH outDataset;
    const char* outname = fname.c_str();

    int nrows, ncols;
    int row, col;
    ncols = src.xsize;
    nrows = src.ysize;

    int count;
    count = 0;
    while (count < DIMGT) {
        geotransform[count] = src.geotransform[count];
        count++;
    }
    outHandleDriver = GDALGetDriverByName("GTiff");
    outDataset = GDALCreate(outHandleDriver, outname, ncols, nrows, 4, GDT_UInt16, NULL);
    GDALSetGeoTransform(outDataset, geotransform);
    GDALSetProjection(outDataset, projection);

    uint16_t* scanLineblue = (uint16_t*)CPLMalloc(sizeof(uint16_t) * ncols);
    uint16_t* scanLinegreen = (uint16_t*)CPLMalloc(sizeof(uint16_t) * ncols);
    uint16_t* scanLinered = (uint16_t*)CPLMalloc(sizeof(uint16_t) * ncols);
    uint16_t* scanLineNIR = (uint16_t*)CPLMalloc(sizeof(uint16_t) * ncols);

    GDALRasterBandH handleRedBand;
    GDALRasterBandH handleGreenBand;
    GDALRasterBandH handleBlueBand;
    GDALRasterBandH handleNIRBand;

    handleRedBand = GDALGetRasterBand(outDataset, 1);
    handleGreenBand = GDALGetRasterBand(outDataset, 2);
    handleBlueBand = GDALGetRasterBand(outDataset, 3);
    handleNIRBand = GDALGetRasterBand(outDataset, 4);

    for (row = 0; row < nrows; row++) {
        for (col = 0; col < ncols; col++) {
            scanLinered[col] = rband[row][col];
            scanLinegreen[col] = gband[row][col];
            scanLineblue[col] = bband[row][col];
            scanLineNIR[col] = nband[row][col];

        }
        GDALRasterIO(handleRedBand, GF_Write, 0, row, ncols, 1, scanLinered, ncols, 1, GDT_UInt16, 0, 0);
        GDALRasterIO(handleGreenBand, GF_Write, 0, row, ncols, 1, scanLinegreen, ncols, 1, GDT_UInt16, 0, 0);
        GDALRasterIO(handleBlueBand, GF_Write, 0, row, ncols, 1, scanLineblue, ncols, 1, GDT_UInt16, 0, 0);
        GDALRasterIO(handleNIRBand, GF_Write, 0, row, ncols, 1, scanLineNIR, ncols, 1, GDT_UInt16, 0, 0);
    }
    GDALClose(outDataset);
    return 0;
}

int writeGeotiff(
    Geotiff pan,
    float** red,
    float** green,
    float** blue,
    float** NIR,
    string& fname
) {

    char* projection = pan.projection;
    double geotransform[DIMGT];
    GDALDriverH outHandleDriver;
    GDALDatasetH outDataset;
    int nrows, ncols;
    int row, col;
    ncols = pan.xsize;
    nrows = pan.ysize;
    const char* outname = fname.c_str();



    int count;
    count = 0;
    while (count < DIMGT) {
        geotransform[count] = pan.geotransform[count];
        count++;
    }

    outHandleDriver = GDALGetDriverByName("GTiff");
    outDataset = GDALCreate(outHandleDriver,
        outname,
        ncols, nrows, 4,
        GDT_Float32, NULL);
    GDALSetGeoTransform(outDataset, geotransform);
    GDALSetProjection(outDataset, projection);


    float* scanLineRed = (float*)CPLMalloc(sizeof(float) * ncols);
    float* scanLineGreen = (float*)CPLMalloc(sizeof(float) * ncols);
    float* scanLineBlue = (float*)CPLMalloc(sizeof(float) * ncols);
    float* scanLineNIR = (float*)CPLMalloc(sizeof(float) * ncols);


    GDALRasterBandH handleRedBand;
    GDALRasterBandH handleGreenBand;
    GDALRasterBandH handleBlueBand;
    GDALRasterBandH handleNIRBand;

    handleRedBand = GDALGetRasterBand(outDataset, 1);
    handleGreenBand = GDALGetRasterBand(outDataset, 2);
    handleBlueBand = GDALGetRasterBand(outDataset, 3);
    handleNIRBand = GDALGetRasterBand(outDataset, 4);

    for (row = 0; row < nrows; row++) {

        for (col = 0; col < ncols; col++) {
            scanLineRed[col] = red[row][col];
            scanLineGreen[col] = green[row][col];
            scanLineBlue[col] = blue[row][col];
            scanLineNIR[col] = NIR[row][col];
        }

        GDALRasterIO(
            handleRedBand, GF_Write, 0, row, ncols, 1,
            scanLineRed, ncols, 1,
            GDT_Float32, 0, 0
        );

        GDALRasterIO(
            handleGreenBand, GF_Write, 0, row, ncols, 1,
            scanLineGreen, ncols, 1,
            GDT_Float32, 0, 0
        );

        GDALRasterIO(
            handleBlueBand, GF_Write, 0, row, ncols, 1,
            scanLineBlue, ncols, 1,
            GDT_Float32, 0, 0
        );

        GDALRasterIO(
            handleNIRBand, GF_Write, 0, row, ncols, 1,
            scanLineNIR, ncols, 1,
            GDT_Float32, 0, 0
        );


    }

    GDALClose(outDataset);
    return 0;
}


int TOA_reflectance(const char* img_bands[], string& metadata, string& result_path) {

    Geotiff r_band, g_band, b_band, nir_band;

    b_band = readGeotiff(img_bands[0]);
    g_band = readGeotiff(img_bands[1]);
    r_band = readGeotiff(img_bands[2]);
    nir_band = readGeotiff(img_bands[3]);

    int ncol = b_band.xsize;
    int nrow = b_band.ysize;

    float** b_radi = b_band.band;
    float** g_radi = g_band.band;
    float** r_radi = r_band.band;
    float** nir_radi = nir_band.band;


    GeotiffMeta b_meta, g_meta, r_meta, n_meta;
    for (int i = 0; i < 4; i++) {
        if (i == 0) {

            b_meta = readGeoTiffMeta(metadata, i + 2);

        }
        else if (i == 1) {

            g_meta = readGeoTiffMeta(metadata, i + 2);

        }
        else if (i == 2) {

            r_meta = readGeoTiffMeta(metadata, i + 2);

        }
        else if (i == 3) {

            n_meta = readGeoTiffMeta(metadata, i + 2);

        }
    }

    float** r_TOA_refl = (float**)malloc(nrow * sizeof(float*));
    float** g_TOA_refl = (float**)malloc(nrow * sizeof(float*));
    float** b_TOA_refl = (float**)malloc(nrow * sizeof(float*));
    float** n_TOA_refl = (float**)malloc(nrow * sizeof(float*));

    for (int row = 0; row < nrow; row++) {
        r_TOA_refl[row] = (float*)malloc(ncol * sizeof(float));
        g_TOA_refl[row] = (float*)malloc(ncol * sizeof(float));
        b_TOA_refl[row] = (float*)malloc(ncol * sizeof(float));
        n_TOA_refl[row] = (float*)malloc(ncol * sizeof(float));

        for (int col = 0; col < ncol; col++) {
            r_TOA_refl[row][col] = ((r_radi[row][col] * r_meta.refl_mul) + r_meta.refl_add);
            g_TOA_refl[row][col] = ((g_radi[row][col] * g_meta.refl_mul) + g_meta.refl_add);
            b_TOA_refl[row][col] = ((b_radi[row][col] * b_meta.refl_mul) + b_meta.refl_add);
            n_TOA_refl[row][col] = ((nir_radi[row][col] * n_meta.refl_mul) + n_meta.refl_add);


        }
    }
    writeGeotiff(b_band, r_TOA_refl, g_TOA_refl, b_TOA_refl, n_TOA_refl, result_path);
    return 0;

}