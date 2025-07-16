from gdal_class import gdal_data
import numpy as np
from osgeo import gdal
import sys


def resample_img(dataset, GSD):
    x_gsd = GSD[0]
    y_gsd = GSD[1]

    resampled_ds = gdal.Warp(
        '',
        dataset,
        format='MEM',
        xRes=x_gsd,
        yRes=y_gsd,
        resampleAlg='cubic'
    )
    return resampled_ds


def run_resample(src_path, p_path):
    p_data = gdal.Open(p_path)
    p_geoinfo = (p_data.GetGeoTransform()[1], p_data.GetGeoTransform()[5])

    tmp_dataset = gdal.Open(src_path)
    tmp_resampled = resample_img(tmp_dataset, p_geoinfo)
    tmp_gdal_data = gdal_data.from_dataset(tmp_resampled)

    last_slash = max(src_path.rfind('/'), src_path.rfind('\\'))
    tmp_gdal_data.filename = src_path[last_slash + 1:]

    return tmp_gdal_data

### main code
src_path = str(sys.argv[1])#기존 해상도
dst_path =  str(sys.argv[2])#최종 해상도
result_path = r'D:\project\RS_img_preprocess\data\output\resampling'

result_dataset = run_resample(src_path, dst_path)
result_dataset.Write(result_path, 'resampled')













