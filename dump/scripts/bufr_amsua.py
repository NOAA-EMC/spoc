#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions

AMUA_1B_MAPPING = map_path('bufr_amsua_1bamua.yaml')
AMUA_ES_MAPPING = map_path('bufr_amsua_esamua.yaml')


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
        self.a_sp = self.a_space * COSMIC_BACKGROUND_TEMP


class PressureObsBuilder(ObsBuilder):
    def __init__(self):
        map_dict = {'1bamua': AMUA_1B_MAPPING,
                    'esamua': AMUA_ES_MAPPING}

        super().__init__(map_dict, log_name=os.path.basename(__file__))

    def _remove_ant_corr(self, i, ac, ifov, t):
        # AC:             Structure containing the antenna correction coefficients for the sensor of interest.
        # iFOV:           The FOV index for a scanline of the sensor of interest.
        # T:              On input, this argument contains the brightness

        t = ac.a_ep[i, ifov] * t + ac.a_sp[i, ifov]
        t[(ifov < 1) | (ifov > ac.n_fovs)] = [INVALID]
        return t

    def _apply_ant_corr(self, i, ac, ifov, t):
        # t:              on input, this argument contains the antenna temperatures for the sensor channels.
        t = (t - ac.a_sp[i, ifov]) / ac.a_ep[i, ifov]
        t[(ifov < 1) | (ifov > ac.n_fovs)] = [INVALID]
        return t

    def _apply_corr(self, sat_id, ta, ifov):
        ac = ACCoeff(nc_dir, sat_id=sat_id)
        if sat_id not in ['n15', 'n16']:
            # Convert antenna temperature to brightness temperature
            ifov = ifov.astype(int) - 1
            for i in range(ta.shape[1]):
                self.log.debug(f'inside loop for allpy ta to tb: i = {i}')
                x = ta[:, i]
                if YAML_NORMAL:
                    x = _apply_ant_corr(i, ac, ifov, x)
                else:
                    x = _remove_ant_corr(i, ac, ifov, x)
                x[x >= INVALID] = INVALID
                ta[:, i] = x
        return ta

    def _re_map_variable(self, comm, container):
        # read_bufrtovs.f90
        # antcorr_application.f90
        # search the keyword “ta2tb” for details
        sat_ids = container.all_sub_categories()
        for sat_id in sat_ids:
            self.log.info(f'Converting for {sat_id[0]}, ...')
            ta = container.get('variables/brightnessTemperature', sat_id)
            if ta.shape[0]:
                ifov = container.get('variables/fieldOfViewNumber', sat_id)
                tb = _apply_corr(comm, sat_id[0], ta, ifov)
                container.replace('variables/brightnessTemperature', tb, sat_id)

    def make_obs(self, comm, input_path):
        cache = bufr.DataCache.has(input_path, yaml_path)
        if cache:
            logging(comm, f'The cache existed get data container from it')
            container = bufr.DataCache.get(input_path, yaml_path)
        else:
            # If cacache does not exist, get data into cache
            # Get data info container first
            logging(comm, f'The cache is not existed')
            container = bufr.Parser(input_path, yaml_path).parse()
        return cache, container


def create_obs_group(input_path1, input_path2, yaml_1b, yaml_es, category, env):
    comm = bufr.mpi.Comm(env["comm_name"])

    logging(comm, f'Imput_path: {input_path1}, {input_path2}, and category: {category}')
    logging(comm, f'Entering function to create obs group for {category} with yaml path {yaml_es} and {yaml_1b}')
    cache_1, container_1 = _make_obs(comm, input_path1, yaml_es)
    cache_2, container_2 = _make_obs(comm, input_path2, yaml_1b)

    container = container_1

    if not cache_2:
        logging(comm, f'No chahe, remap and append it')
        _re_map_variable(comm, container_2)
        container.append(container_2)
        logging(comm, 'Container append done')

    data = Encoder(_make_description(yaml_es)).encode(container)[(category,)]
    _mark_one_data(comm, cache_1, input_path1, yaml_es, category, container=container_1)
    _mark_one_data(comm, cache_2, input_path2, yaml_1b, category, container=container_2)
    return data


def create_obs_file(input_path1, input_path2, yaml_1b, yaml_es, output_path):

    comm = bufr.mpi.Comm("world")

    logging(comm, f'Imput_path: {input_path1} and {input_path2}')
    logging(comm, f'Entering function with yaml path {yaml_es} and {yaml_1b}')
    cache_1, container_1 = _make_obs(comm, input_path1, yaml_es)
    cache_2, container_2 = _make_obs(comm, input_path2, yaml_1b)

    _re_map_variable(comm, container_2)

    container = container_1
    container.append(container_2)
    logging(comm, 'Container append done')

    description = _make_description(yaml_es)

    # Encode the data
    if comm.rank() == 0:
        netcdfEncoder(description).encode(container, output_path)

    logging(comm, f'Return the encoded data')


if __name__ == '__main__':
    start_time = time.time()
    bufr.mpi.App(sys.argv)
    comm = bufr.mpi.Comm("world")

    # Required input arguments as positional arguments
    parser = argparse.ArgumentParser(description="Convert BUFR to NetCDF using a mapping file.")
    parser.add_argument('input_path1', type=str, help='Input BUFR file for 1b')
    parser.add_argument('input_path2', type=str, help='Input BUFR file for es')
    parser.add_argument('yaml_1b', type=str, help='BUFR2IODA Mapping File for 1b')
    parser.add_argument('yaml_es', type=str, help='BUFR2IODA Mapping File for es')
    parser.add_argument('output', type=str, help='Output NetCDF file')

    args = parser.parse_args()
    input_path1 = args.input_path1
    input_path2 = args.input_path2
    yaml_1b = args.yaml_1b
    yaml_es = args.yaml_es
    output = args.output

    create_obs_file(input_path1, input_path2, yaml_1b, yaml_es, output)

    end_time = time.time()
    running_time = end_time - start_time
    logging(comm, f'Total running time: {running_time}')
