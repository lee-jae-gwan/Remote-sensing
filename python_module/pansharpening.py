from gdal_class import gdal_data
from osgeo import gdal
import numpy as np
import glob
from sklearn.linear_model import LinearRegression

def resample_img(dataset, GSD, outbounds):
    x_gsd = GSD[0]
    y_gsd = GSD[1]

    resampled_ds = gdal.Warp(
        '',
        dataset,
        format='MEM',
        xRes=x_gsd,
        yRes=y_gsd,
        outputBounds = outbounds,
        resampleAlg='cubic'
    )
    return resampled_ds


def run_resample(src_path, p_path):
    p_data = gdal.Open(p_path)
    p_geoinfo = (p_data.GetGeoTransform()[1], p_data.GetGeoTransform()[5])
    outbounds = p_data.GetGeoTransform()

    x_min = outbounds[0]
    y_max = outbounds[3]
    x_max = x_min + p_data.RasterXSize * outbounds[1]
    y_min = y_max + p_data.RasterYSize * outbounds[5]

    outbounds = (x_min, y_min, x_max, y_max)

    tmp_dataset = gdal.Open(src_path)
    tmp_resampled = resample_img(tmp_dataset, p_geoinfo, outbounds)
    tmp_gdal_data = gdal_data.from_dataset(tmp_resampled)

    last_slash = max(src_path.rfind('/'), src_path.rfind('\\'))
    tmp_gdal_data.filename = src_path[last_slash + 1:]

    return tmp_gdal_data

def estimate_coefficients(ms_low, pan_low):
    H, W, C = ms_low.shape
    X = ms_low.reshape(-1, C)  # (H*W, 4)
    y = pan_low.reshape(-1, 1)  # (H*W, 1)
    reg = LinearRegression().fit(X, y)
    return reg.coef_[0].astype('float32'), reg.intercept_.astype('float32')


def make_GIHSA_pansharpen_img(paths, resampled_datas):
    hir_pan = gdal_data.Open(paths['pan'])
    hir_pan_mask = hir_pan.band > 0

    low_pan = run_resample(paths['pan'], paths['ms'][0])
    low_pan = low_pan.band

    ms_bands = np.dstack([
        gdal_data.Open(path).band
        for path in paths['ms']
    ])

    resampled_ms_bands = list(
        [resampled_datas[0].band, resampled_datas[1].band, resampled_datas[2].band, resampled_datas[3].band])
    nrow, ncol = resampled_datas[0].band.shape
    hir_pan.band = hir_pan.band[:nrow, :ncol]
    hir_pan_mask = hir_pan_mask[:nrow, :ncol]

    resam_data = [
        np.where(hir_pan_mask > 0, band.band, 0)
        for band in resampled_datas
    ]

    w, b = estimate_coefficients(ms_bands, low_pan)
    intensity = sum(w[i] * resampled_datas[i].band for i in range(len(resampled_datas))) + b

    delta = np.where(hir_pan_mask, hir_pan.band - intensity, 0)

    fused = [band.band + delta for band in resampled_datas]

    return fused


def write_bands(datasets, dst_path, task):
    for i in range(len(datasets)):
        datasets[i].Write(dst_path, task)


def run_pansharpening(band_paths, p_path, result_path):
    data_path = {'ms': [], 'pan': []}
    for i in band_paths.values():
        data_path['ms'].append(i)

    data_path['pan'] = p_path

    pan_dataset = gdal_data.Open(data_path['pan'])

    resampled_ms_dataset = []
    print("----------------------start resampleing-------------------------------------")
    for i in range(len(data_path['ms'])):
        resampled_ms_dataset.append(run_resample(data_path['ms'][i], data_path['pan']))

    print("----------------------start GIHS-------------------------------------")
    GIHSA_bands = make_GIHSA_pansharpen_img(data_path, resampled_ms_dataset)

    for i in range(len(resampled_ms_dataset)):
        tmp_band = GIHSA_bands[i]
        min_val = GIHSA_bands[i].min()
        if min_val < 0:
            resampled_ms_dataset[i].band = (tmp_band - min_val).astype('uint16')
        else:
            resampled_ms_dataset[i].band = (tmp_band + min_val).astype('uint16')
        resampled_ms_dataset[i].Write(result_path, '')