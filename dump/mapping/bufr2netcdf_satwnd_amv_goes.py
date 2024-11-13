import sys 
import os 
import argparse
import numpy as np
import bufr
import time
import warnings
from bufr.encoders import netcdf
from wxflow import Logger

def get_description(yaml_path, update=False):

    description = bufr.encoders.Description(yaml_path)

    if update:

        description.add_variable(name='ObsType/windEastward',
                                 source='variables/obstype_uwind',
                                 units='1',
                                 longName='Observation Type based on Satellite-derived Wind Computation Method and Spectral Band')

        description.add_variable(name='ObsType/windNorthward',
                                 source='variables/obstype_vwind',
                                 units='1',
                                 longName='Observation Type based on Satellite-derived Wind Computation Method and Spectral Band')

        description.add_variable(name='ObsValue/windEastward',
                                 source='variables/windEastward',
                                 units='m/s',
                                 longName='Eastward Wind Component')

        description.add_variable(name='ObsValue/windNorthward',
                                 source='variables/windNorthward',
                                 units='m/s',
                                 longName='Northward Wind Component')

    return description

def get_all_keys(nested_dict):

    keys = []
    for key, value in nested_dict.items():
        keys.append(key)
        if isinstance(value, dict):
            keys.extend(get_all_keys(value))
    return keys

def Compute_WindComponents_from_WindDirection_and_WindSpeed(wdir, wspd):

    uob = (-wspd * np.sin(np.radians(wdir))).astype(np.float32)
    vob = (-wspd * np.cos(np.radians(wdir))).astype(np.float32)

    return uob, vob

def Get_ObsType(swcm, chanfreq):

    obstype = swcm.copy()

    # Use numpy vectorized operations
    obstype = np.where(swcm == 5, 247, obstype)  # WVCA/DL
    obstype = np.where(swcm == 3, 246, obstype)  # WVCT
    obstype = np.where(swcm == 2, 251, obstype)  # VIS
    obstype = np.where(swcm == 1, 245, obstype)  # IRLW

    condition = np.logical_and(swcm == 1, chanfreq >= 50000000000000.0)  # IRSW
    obstype = np.where(condition, 240, obstype)

    if not np.any(np.isin(obstype, [247, 246, 251, 245, 240])):
        raise ValueError("Error: Unassigned ObsType found ... ")

    return obstype

def create_obs_group(yaml_path, input_path, output_path, logger):

    # Get MPI communicator
    comm = bufr.mpi.Comm("world")
    rank = comm.rank()

    # Generate keys for cache
    inputKey = input_path
    mappingKey = yaml_path

    if rank == 0: logger.info('Input parameters')
    if rank == 0: logger.info(f' ... input file: {input_path}')
    if rank == 0: logger.info(f' ... output file: {output_path}')
    if rank == 0: logger.info(f' ... mapping file: {yaml_path}')

    # Get container from mapping file first
    if rank == 0: logger.info('Create data container based on mapping file')
    container = bufr.Parser(input_path, yaml_path).parse(comm)
    if rank == 0: logger.info(f' ... all_sub_categories =  {container.all_sub_categories()}')
    if rank == 0: logger.debug(' ... container list (original):\n' + '\n'.join(str(item) for item in container.list()))

    categories = container.all_sub_categories()
    if rank == 0: logger.info('Loop through data categories to add new/derived variables to data container')
    for cat in categories:  

        description = get_description(yaml_path, update=True)
        if rank == 0: logger.info(f' ... category = {cat[0]}')

        satid = container.get('variables/satelliteId', cat)
        if satid.size == 0:
            if rank == 0: logger.info(f' ... category {cat[0]} does not exist in input file')

            paths = container.get_paths('variables/windComputationMethod', cat)
            obstype = container.get('variables/windComputationMethod', cat)
            container.add('variables/obstype_uwind', obstype, paths, cat)
            container.add('variables/obstype_vwind', obstype, paths, cat)

            paths = container.get_paths('variables/windSpeed', cat)
            wob = container.get('variables/windSpeed', cat)
            container.add('variables/windEastward', wob, paths, cat)
            container.add('variables/windNorthward', wob, paths, cat)

        else:
            # Add new variables: ObsType/windEastward & ObsType/windNorthward 
            swcm = container.get('variables/windComputationMethod', cat)
            chanfreq = container.get('variables/sensorCentralFrequency', cat)

            if rank == 0: logger.debug(f' ... swcm min/max = {swcm.min()} {swcm.max()}')
            if rank == 0: logger.debug(f' ... chanfreq min/max =  {chanfreq.min()} {chanfreq.max()}')

            obstype = Get_ObsType(swcm, chanfreq)

            if rank == 0: logger.info(f' ... adding obstype to data container')
            if rank == 0: logger.debug(f' ... obstype = {obstype}')
            if rank == 0: logger.debug(f' ... obstype min/max =  {obstype.min()} {obstype.max()}')

            paths = container.get_paths('variables/windComputationMethod', cat)
            container.add('variables/obstype_uwind', obstype, paths, cat)
            container.add('variables/obstype_vwind', obstype, paths, cat)

            # Add new variables: ObsValue/windEastward & ObsValue/windNorthward 
            wdir = container.get('variables/windDirection', cat)
            wspd = container.get('variables/windSpeed', cat)

            if rank == 0: logger.debug(f' ... wdir min/max = {wdir.min()} {wdir.max()}')
            if rank == 0: logger.debug(f' ... wspd min/max = {wspd.min()} {wspd.max()}')

            uob, vob = Compute_WindComponents_from_WindDirection_and_WindSpeed(wdir, wspd)

            if rank == 0: logger.info(f' ... adding windEastward and windNorthward to data container')
            if rank == 0: logger.debug(f' ... uob min/max = {uob.min()} {uob.max()}')
            if rank == 0: logger.debug(f' ... vob min/max = {vob.min()} {vob.max()}')

            paths = container.get_paths('variables/windSpeed', cat)
            container.add('variables/windEastward', uob, paths, cat)
            container.add('variables/windNorthward', vob, paths, cat)

    # Check
    if rank == 0: logger.debug(' container list (updated):\n' + '\n'.join(str(item) for item in container.list()))

    # Gather the DataContainer data from all MPI ranks [Optional]
    if rank == 0: logger.info('Gather data container from all MPI tasks rank 0 task')
    container.gather(comm)   

    # Encode all categories
    if rank == 0:
        logger.info('Encoding the data at MPI rank 0')
        netcdf.Encoder(description).encode(container, output_path)


if __name__ == '__main__':

    start_time = time.time()

    bufr.mpi.App(sys.argv) 
    comm = bufr.mpi.Comm("world")
    rank = comm.rank()           

    # Required input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input NetCDF', required=True)
    parser.add_argument('-m', '--mapping', type=str, help='Mapping File', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output NetCDF', required=True)
    parser.add_argument('-v', '--verbose', help='print debug logging info', action='store_true')

    args = parser.parse_args()
    mapping = args.mapping
    infile = args.input
    output = args.output

    log_level = 'DEBUG' if args.verbose else 'INFO'
    logger = Logger('BUFR2IODA_satwind_amv_goes.py', level=log_level, colored_log=True)    

    create_obs_group(mapping, infile, output, logger)

    end_time = time.time()
    running_time = end_time - start_time
    if rank == 0: logger.info(f'Total running time: {running_time}')
