from osgeo import gdal
import numpy as np
from sklearn.linear_model import LinearRegression
from gdal_class import gdal_data
import glob

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


def calc_img_size(datasets):
    smallest = min(datasets, key=lambda x: x.band.shape[0] * x.band.shape[1])
    return smallest.band.shape


def match_img_size(datasets, rows, cols):
    if isinstance(datasets, list):
        for item in datasets:
            item.band = item.band[:rows, :cols]
    elif isinstance(datasets, gdal_data):
        datasets.band = datasets.band[:rows, :cols]
    else:
        raise TypeError("지원하지 않는 데이터 타입입니다.")

    return datasets


def create_valid_band(datasets):
    tmp_min = float('inf')
    mask = np.zeros_like(datasets[0].band)

    for i in range(len(datasets)):
        tmp_band = datasets[i].band
        nopad_pixels = len(tmp_band[tmp_band > 0])

        if tmp_min > nopad_pixels:
            tmp_min = nopad_pixels
            mask = tmp_band > 0

    b = datasets[0].band[mask]
    g = datasets[1].band[mask]
    r = datasets[2].band[mask]
    n = datasets[3].band[mask]
    p = datasets[4].band[mask]

    b[b == 0] = 1
    g[g == 0] = 1
    r[r == 0] = 1
    n[n == 0] = 1
    p[p == 0] = 1

    return b, g, r, n, p, mask


def brovey_pansharpening(b, g, r, n, p):
    sum_band = b + g + r + n

    sum_band[sum_band == 0] = 1e-8

    pan_b = (b / sum_band) * p
    pan_g = (g / sum_band) * p
    pan_r = (r / sum_band) * p
    pan_n = (n / sum_band) * p
    return pan_b, pan_g, pan_r, pan_n


def make_dataset(path):
    dataset = {'ms': [], 'pan': []}
    file_list = glob.glob(path + '/*.tif')

    for file in file_list:
        if '_PS' not in file and '_P' not in file:
            dataset['ms'].append(file)
        if '_P.' in file:
            dataset['pan'].append(file)
    return dataset


def make_brovey_pansharpen_img(dataset, pan_img):
    min_row, min_col = calc_img_size(dataset)
    match_datasets = match_img_size(dataset, min_row, min_col)
    pan_datasets = match_img_size(pan_img, min_row, min_col)

    match_datasets.append(pan_datasets)

    valid_b, valid_g, valid_r, valid_n, valid_p, mask = create_valid_band(match_datasets)

    pan_b, pan_g, pan_r, pan_n = brovey_pansharpening(valid_b, valid_g, valid_r, valid_n, valid_p)

    pan_bands = [pan_b, pan_g, pan_r, pan_n]
    for i, pan in enumerate(pan_bands):
        match_datasets[i].band[mask] = pan

    return match_datasets


def estimate_coefficients(ms_low, pan_low):
    H, W, C = ms_low.shape
    X = ms_low.reshape(-1, C)  # (H*W, 4)
    y = pan_low.reshape(-1, 1)  # (H*W, 1)
    reg = LinearRegression().fit(X, y)
    return reg.coef_[0].astype('float32'), reg.intercept_.astype('float32')


def make_GIHSA_pansharpen_img(p_data, resampled_datas, low_path):
    p_data_mask = p_data.band > 0

    low_pan = run_resample(p_data['pan'][0], low_path)
    low_pan = low_pan.band

    resampled_ms_bands = list([resampled_datas['R'].band, resampled_datas['G'].band, resampled_datas['B'].band, resampled_datas['N'].band])
    nrow, ncol = resampled_datas[0].band.shape
    p_data.band = p_data.band[:nrow, :ncol]
    p_data_mask = p_data_mask[:nrow, :ncol]

    resam_data = [np.where(p_data_mask > 0, band.band, 0) for band in resampled_datas]

    w, b = estimate_coefficients(resampled_ms_bands, low_pan)
    intensity = sum(w[i] * resampled_datas[i].band for i in range(len(resampled_datas))) + b

    delta = np.where(hir_pan_mask, hir_pan.band - intensity, 0)

    gihsa_pan = [band.band + delta for band in resampled_datas]

    return gihsa_pan


def write_bands(datasets, dst_path, task):
    for i in range(len(datasets)):
        datasets[i].Write(dst_path, task)
