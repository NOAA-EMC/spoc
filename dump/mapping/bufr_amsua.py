#!/usr/bin/env python3
import bufr
import netCDF4 as nc
import os
from pyioda.ioda.Engines.Bufr import Encoder
from wxflow import Logger

logger = Logger(os.path.basename(__file__), level='INFO', colored_log=True)

YAML_NORMAL = False  # current as normal need remap for path2/1bama

R1000 = 1000.0
R1000000 = 1000000.0
INVALID = R1000
# Cosmic background temperature. Taken from Mather,J.C. et. al., 1999, "Calibrator Design for the COBE
# Far-Infrared Absolute Spectrophotometer (FIRAS)"Astrophysical Journal, vol 512, pp 511-520
TSPACE = 2.7253

nc_dir = './aux'


class ACCoeff:
    def __init__(self, ac_dir, sat_id='n19'):
        file_name = os.path.join(ac_dir, 'amsua_' + sat_id + '.ACCoeff.nc')
        nc_file = nc.Dataset(file_name)
        self.n_fovs = len(nc_file.dimensions['n_FOVs'])
        self.n_channels = len(nc_file.dimensions['n_Channels'])
        self.a_earth = nc_file.variables['A_earth'][:]
        self.a_platform = nc_file.variables['A_platform'][:]
        self.a_space = nc_file.variables['A_space'][:]
        self.a_ep = self.a_earth + self.a_platform
        self.a_sp = self.a_space * TSPACE


def remove_ant_corr(i, ac, ifov, t):
    # AC:             Structure containing the antenna correction coefficients for the sensor of interest.
    # iFOV:           The FOV index for a scanline of the sensor of interest.
    # T:              On input, this argument contains the brightness

    t = ac.a_ep[i, ifov] * t + ac.a_sp[i, ifov]
    t[(ifov < 1) | (ifov > ac.n_fovs)] = [INVALID]
    return t


def apply_ant_corr(i, ac, ifov, t):
    # t:              on input, this argument contains the antenna temperatures for the sensor channels.
    t = (t - ac.a_sp[i, ifov]) / ac.a_ep[i, ifov]
    t[(ifov < 1) | (ifov > ac.n_fovs)] = [INVALID]
    return t


def get_description(yaml_path):
    description = bufr.encoders.Description(yaml_path)
    return description


def apply_corr(sat_id, ta, ifov):
    ac = ACCoeff(nc_dir, sat_id=sat_id)
    if sat_id not in ['n15', 'n16']:
        # Convert antenna temperature to brightness temperature
        ifov = ifov.astype(int) - 1
        for i in range(ta.shape[1]):
            logger.debug(f'inside loop for allpy ta to tb: i = {i}')
            x = ta[:, i]
            if YAML_NORMAL:
                x = apply_ant_corr(i, ac, ifov, x)
            else:
                x = remove_ant_corr(i, ac, ifov, x)
            x[x >= R1000] = R1000000
            ta[:, i] = x
    return ta


def re_map_variable(container):
    # read_bufrtovs.f90
    # antcorr_application.f90
    # search the keyword “ta2tb” for details
    sat_ids = container.all_sub_categories()
    for sat_id in sat_ids:
        logger.info(f'Converting for {sat_id[0]}, ...')
        ta = container.get('variables/brightnessTemperature', sat_id)
        if ta.shape[0]:
            ifov = container.get('variables/fieldOfViewNumber', sat_id)
            tb = apply_corr(sat_id[0], ta, ifov)
            container.replace('variables/brightnessTemperature', tb, sat_id)


def get_one_data(input_path, yaml_path, category):
    cache = bufr.DataCache.has(input_path, yaml_path)
    if cache:
        logger.info(f'The cache existed get data container from it')
        container = bufr.DataCache.get(input_path, yaml_path)
    else:
        # If cacache does not exist, get data into cache
        # Get data info container first
        logger.info(f'The cache is not existed')
        container = bufr.Parser(input_path, yaml_path).parse()
    return cache, container


def mark_one_data(cache, input_path, yaml_path, category, container=None):
    if cache:
        logger.info(f'The cache existed get data container from it')
        bufr.DataCache.mark_finished(input_path, yaml_path, [category])
    else:
        logger.info(f'add original container list into a cache = {container.list()}')
        bufr.DataCache.add(input_path, yaml_path, container.all_sub_categories(), container)
        bufr.DataCache.mark_finished(input_path, yaml_path, [category])


def create_obs_group(input_path1, input_path2, yaml_1b, yaml_es, category):
    logger.info(f'imput_path: {input_path1}, {input_path2}, and category: {category}')
    logger.info(f'Entering function to create obs group for {category} with yaml path {yaml_es} and {yaml_1b}')
    cache_1, container_1 = get_one_data(input_path1, yaml_es, category)
    cache_2, container_2 = get_one_data(input_path2, yaml_1b, category)

    re_map_variable(container_2)

    container = container_1
    container.append(container_2)

    logger.info('Container append done')
    data = Encoder(get_description(yaml_es)).encode(container)[(category,)]
    mark_one_data(cache_1, input_path1, yaml_es, category, container=container_1)
    mark_one_data(cache_2, input_path2, yaml_1b, category, container=container_2)
    return data

