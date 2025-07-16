from gdal_class import gdal_data
from datetime import datetime
import xml.etree.ElementTree as ET
import sys
import ee
import math
from Py6S import*
from google.oauth2 import service_account

from atmospheric import Atmospheric


def spectralResponseFunction(bandname):
    """
    Extract spectral response function for given band name
    """
    bandSelect = {
        'B1': PredefinedWavelengths.S2A_MSI_01,
        'B2': PredefinedWavelengths.S2A_MSI_02,
        'B3': PredefinedWavelengths.S2A_MSI_03,
        'B4': PredefinedWavelengths.S2A_MSI_04,
        'B5': PredefinedWavelengths.S2A_MSI_05,
        'B6': PredefinedWavelengths.S2A_MSI_06,
        'B7': PredefinedWavelengths.S2A_MSI_07,
        'B8': PredefinedWavelengths.S2A_MSI_08,
        'B8A': PredefinedWavelengths.S2A_MSI_8A,
        'B9': PredefinedWavelengths.S2A_MSI_09,
        'B10': PredefinedWavelengths.S2A_MSI_10,
        'B11': PredefinedWavelengths.S2A_MSI_11,
        'B12': PredefinedWavelengths.S2A_MSI_12,
    }
    return Wavelength(bandSelect[bandname])


def get_date(meta_data):
    for tmp_date in meta_data.iter('ImagingStartTime'):
        date = tmp_date.find('UTC').text

    Y = date[:4]
    m = date[4:6]
    d = date[6:8]
    H = date[8:10]
    M = date[10:12]
    S = date[12:14]
    ee_date = ee.Date(Y + '-' + m + '-' + d)

    dt_part = date[:14]
    frac_part = date[15:]

    microsecond = int(frac_part[:6].ljust(6, '0'))
    dt = datetime.strptime(dt_part, '%Y%m%d%H%M%S').replace(microsecond=microsecond)

    return ee_date, dt


def get_geometry(meta_data):
    for child in meta_data.iter('ImageGeogCenter'):
        latitude = float(child.find('Latitude').text)
        longitude = float(child.find('Longitude').text)

    return latitude, longitude


def get_atmos(geom, date):
    h2o = Atmospheric.water(geom, date).getInfo()
    o3 = Atmospheric.ozone(geom, date).getInfo()
    aot = Atmospheric.aerosol(geom, date).getInfo()

    return h2o, o3, aot


def get_solar_zenith(meta_data):
    sun_ang = next(meta_data.iter("SunAngle"))
    solar_zenith = 90 - float(sun_ang.find("Elevation").text)
    return solar_zenith


def get_bandwith(meta):
    for r in meta.iter("MS3"):
        r_b = r.find('Bandwidth').text
        r_bandwidth = (int(r_b.split('-')[0]) + int(r_b.split('-')[1])) / 2

        for i in r.iter('RadianceConversion'):
            r_gain = float(i.find('Gain').text)
            r_offset = float(i.find('Offset').text)

    return r_gain, r_offset


def get_nadir(meta):
    for sen_ang in meta.iter("Angle"):
        offNadir = float(sen_ang.find("OffNadir").text)

    return offNadir


def get_altitude(meta):
    for alti in meta.iter("SatellitePosition"):
        altitude = round(float(alti.find("Altitude").text))

    return altitude


### 메인코드
def transform_atmos(img, xml_data, band_num):
    band = img
    meta_data = xml_data.getroot()
    print(meta_data)

    ee_date, dt = get_date(meta_data)  # 날짜 추출 및 구글어스엔진의 날짜 타입으로 변환
    print(ee_date, dt)
    longitude, latitude = get_geometry(meta_data)  # 위치 좌표 추출
    geom = ee.Geometry.Point(longitude, latitude)  # 구글어스엔진으로 좌표 저장
    h2o, o3, aot = get_atmos(geom, ee_date)  # 구글어스엔진으로 수증기, 오존, 에어로졸 추출
    solar_zenith = get_solar_zenith(meta_data)  # 태양천정각 계산
    offNadir = get_nadir(meta_data)  # 센서 시야각
    altitude = get_altitude(meta_data)  # 태양 고도
    # sat_alti = get_sat_altitude(meta_data)

    ###대기 모델 제작
    s = SixS()

    s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o, o3)
    s.aero_profile = AeroProfile.Continental
    s.aot550 = aot
    s.geometry = Geometry.User()
    s.geometry.view_z = offNadir  # 모델에 센서 시야각 설정
    s.geometry.solar_z = solar_zenith  # 모델에 태양 천정각 설정
    s.geometry.month = dt.month
    s.geometry.day = dt.day
    s.altitudes.set_sensor_satellite_level()
    # s.altitudes.set_target_custom_altitude(altitude)
    # bandwith = get_bandwith(meta_data)#밴드 대역폭 추출

    s.wavelength = spectralResponseFunction(band_num)  # SRFf를 sentinel 영상에 근사
    s.run()

    Edir = s.outputs.direct_solar_irradiance  # 직접 태양 복사에너지양
    Edif = s.outputs.diffuse_solar_irradiance  # 산란 태양 복사에너지양
    Lp = s.outputs.atmospheric_intrinsic_radiance  # 대기 방출, 산란 복사에너지양
    absorb = s.outputs.trans['global_gas'].upward  # 대기 가스에 의한 상향 투과율 => 가스 흡수 없이 위로 통과할 확률
    scatter = s.outputs.trans['total_scattering'].upward  # 전체 대기 산란 효과를 반영한 => 상향 방향의 총 투과율
    tau2 = absorb * scatter  # 복사에너지가 위로 올라가면서 대기 중 가스 흡수와 산란을 모두 통과할 수 있는 비율

    ref = (band - Lp) * (math.pi) / (tau2 * (Edir + Edif))  # 산출식에 의한 대기보정 결과 저장

    # band.band = ref
    # band.Write(result_path, 'atmos_correction')

    return ref


def run_atmos_correction(img_path, xml_path, json_key_path, task_type, result_path):
    SERVICE_ACCOUNT_KEY = json_key_path
    credentials = service_account.Credentials.from_service_account_file(
        SERVICE_ACCOUNT_KEY,
        scopes=['https://www.googleapis.com/auth/earthengine']
    )

    ee.Initialize(credentials)

    type_bands = {'road': ['B4', 'B3', 'B2','NDWI','B8','NDVI','NDVI_GLCM','NDWI_GLCM'],
                  'bd': ['B4', 'B3', 'B2'],
                  'fire': ['B4', 'B3', 'B2'],
                  'water': ['B3','NDWI','B8'],
                  'land': ['B4', 'B3', 'B2', "B8"]}



    img_dataset = gdal_data.Open(img_path)

    bands = []
    for i in range(img_dataset.band.shape[2]):
        t_band = img_dataset.band[:, :, i]
        bands.append(t_band)

    xml_data = ET.parse(xml_path)

    task_type_band = type_bands[task_type]

    corrected_bands = []
    for i in range(len(task_type_band)):
        if task_type_band[i] in ['NDWI','NDVI','NDVI_GLCM','NDWI_GLCM']:
            corrected_bands.append(bands[i])
        else:
            tmp_corr_band = transform_atmos(bands[i], xml_data, task_type_band[i])
            corrected_bands.append(tmp_corr_band)

    corrected_img = np.dstack(corrected_bands)

    img_dataset.band = corrected_img

    img_dataset.Wrte(result_path, 'AC_img')