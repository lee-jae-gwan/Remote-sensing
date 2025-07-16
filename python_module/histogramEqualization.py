from gdal_class import gdal_data
import numpy as np
import sys
import os

#확률밀도함수 구하기
def calc_pdf(values, max_val):
    pdf = np.zeros(max_val + 1, dtype='int')
    for i in values:
        pdf[i] += 1

    return pdf

#누적분포함수 구하기
def calc_cdf(pdf):
    cdf = np.zeros_like(pdf)
    cdf[0] = pdf[0]

    for i in range(1, (len(pdf))):
        cdf[i] = cdf[i - 1] + pdf[i]

    return cdf

#누적분포함수 정규화=> look-up table 제작
def normalize_cdf(cdf, max_val):
    norm_cdf = np.zeros_like(cdf)

    cdf_min = cdf[cdf > 0].min()
    cdf_max = cdf.max()
    num_pix = cdf[-1]

    for i in range(len(cdf)):
        norm_cdf[i] = round(((cdf[i] - cdf_min) / (num_pix - cdf_min)) * max_val)

    return norm_cdf


def histogramEqualize(dataset):
    src_band = dataset.band.copy()
    equalized_band = dataset.band.copy()

    max_val_dtype = np.iinfo(src_band.dtype).max
    valid_val = src_band[src_band > 0]

    pdf = calc_pdf(valid_val, max_val_dtype + 1)

    cdf = calc_cdf(pdf)

    norm_cdf = normalize_cdf(cdf, max_val_dtype)

    equalized_band = norm_cdf[src_band[src_band > 0]]
    src_band[src_band > 0] = equalized_band.astype(src_band.dtype)

    return src_band

####메인코드
ori_path = sys.argv[1]

histo_path = r'D:\project\RS_img_preprocess\data\output\histogram_equalization'
os.makedirs(histo_path, exist_ok=True)

ori_dataset = gdal_data.Open(ori_path)
equalized_dataset = ori_dataset.copy()
equalized_dataset.band = histogramEqualize(ori_dataset)
equalized_dataset.Write(histo_path, 'histogramEqualized')




