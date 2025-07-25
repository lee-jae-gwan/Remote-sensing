import numpy as np
from gdal_class import gdal_data
import glob

def calc_NDWI(src_path, result_path):

    g_data_path = glob.glob(src_path + '\*_G_pansharpening.tif')
    n_data_path = glob.glob(src_path + '\*_N_pansharpening.tif')

    g_dataset = gdal_data.Open(g_data_path[0])
    n_dataset = gdal_data.Open(n_data_path[0])

    g_band = g_dataset.band.astype('float32')
    n_band = n_dataset.band.astype('float32')


    epsilon = 1e-6 #제로디비전 방지를 위한 매우 작은 수
    ndwi = (g_band - n_band) / (n_band + g_band + epsilon)
    ndwi = np.clip(ndwi, -1, 1)#ndwi 지수를 -1~1 사이의 값으로 지정하기 위한 함수, -1~1사이를 벗어나면 -1과 1로 매핑

    g_mask = g_dataset.band > 0
    ndwi_dataset = g_dataset.copy()
    ndwi_dataset.filename = g_dataset.filename[:g_dataset.filename.index('L2G') + 3] + '_NDWI.tif'
    ndwi_dataset.band = ndwi
    ndwi_dataset.band[~g_mask] = -9999

    ndwi_dataset.Write(result_path, '')





