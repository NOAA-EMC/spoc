#!/usr/bin/env python3
import sys
import os
import argparse
import time
import numpy as np
import bufr
from pyioda.ioda.Engines.Bufr import Encoder as iodaEncoder 
from bufr.encoders.netcdf import Encoder as netcdfEncoder 
from wxflow import Logger

# Initialize Logger
# Get log level from the environment variable, default to 'INFO it not set
log_level = os.getenv('LOG_LEVEL', 'INFO')
logger = Logger('BUFR2IODA_satwnd_amv_goes.py', level=log_level, colored_log=False)

def logging(comm, level, message):
    """
    Logs a message to the console or log file, based on the specified logging level.

    This function ensures that logging is only performed by the root process (`rank 0`) 
    in a distributed computing environment. The function maps the logging level to 
    appropriate logger methods and defaults to the 'INFO' level if an invalid level is provided.

    Parameters:
        comm: object
            The communicator object, typically from a distributed computing framework 
            (e.g., MPI). It must have a `rank()` method to determine the process rank.
        level: str
            The logging level as a string. Supported levels are:
                - 'DEBUG'
                - 'INFO'
                - 'WARNING'
                - 'ERROR'
                - 'CRITICAL'
            If an invalid level is provided, a warning will be logged, and the level 
            will default to 'INFO'.
        message: str
            The message to be logged.

    Behavior:
        - Logs messages only on the root process (`comm.rank() == 0`).
        - Maps the provided logging level to a method of the logger object.
        - Defaults to 'INFO' and logs a warning if an invalid logging level is given.
        - Supports standard logging levels for granular control over log verbosity.

    Example:
        >>> logging(comm, 'DEBUG', 'This is a debug message.')
        >>> logging(comm, 'ERROR', 'An error occurred!')

    Notes:
        - Ensure that a global `logger` object is configured before using this function.
        - The `comm` object should conform to MPI-like conventions (e.g., `rank()` method).
    """

    if comm.rank() == 0:
        # Define a dictionary to map levels to logger methods
        log_methods = {
            'DEBUG': logger.debug,
            'INFO': logger.info,
            'WARNING': logger.warning,
            'ERROR': logger.error,
            'CRITICAL': logger.critical,
        }

        # Get the appropriate logging method, default to 'INFO'
        log_method = log_methods.get(level.upper(), logger.info)

        if log_method == logger.info and level.upper() not in log_methods:
            # Log a warning if the level is invalid
            logger.warning(f'log level = {level}: not a valid level --> set to INFO')

        # Call the logging method
        log_method(message)

def _make_description(mapping_path, update=False):
    description = bufr.encoders.Description(mapping_path)

    if update:
        # Define the variables to be added in a list of dictionaries
        variables = [
            {
                'name': 'ObsType/windEastward',
                'source': 'variables/obstype_uwind',
                'units': '1',
                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
            },
            {
                'name': 'ObsType/windNorthward',
                'source': 'variables/obstype_vwind',
                'units': '1',
                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
            },
            {
                'name': 'ObsValue/windEastward',
                'source': 'variables/windEastward',
                'units': 'm/s',
                'longName': 'Eastward Wind Component',
            },
            {
                'name': 'ObsValue/windNorthward',
                'source': 'variables/windNorthward',
                'units': 'm/s',
                'longName': 'Northward Wind Component',
            },
        ]

        # Loop through each variable and add it to the description
        for var in variables:
            description.add_variable(
                name=var['name'],
                source=var['source'],
                units=var['units'],
                longName=var['longName']
            )

    return description

def compute_wind_components(wdir, wspd):
    """
    Compute the U and V wind components from wind direction and wind speed.

    Parameters:
        wdir (array-like): Wind direction in degrees (meteorological convention: 0° = North, 90° = East).
        wspd (array-like): Wind speed (m/s).

    Returns:
        tuple: U and V wind components as numpy arrays with dtype float32.
    """
    wdir_rad = np.radians(wdir)  # Convert degrees to radians
    u = -wspd * np.sin(wdir_rad)
    v = -wspd * np.cos(wdir_rad)
    
    return u.astype(np.float32), v.astype(np.float32)

def _get_obs_type(swcm, chanfreq):
    """
    Determine the observation type based on `swcm` and `chanfreq`.

    Parameters:
        swcm (array-like): Satellite derived wind calculation method.
        chanfreq (array-like): Satellite channel center frequency (Hz).

    Returns:
        numpy.ndarray: Observation type array.

    Raises:
        ValueError: If any `obstype` is unassigned.
    """

    obstype = swcm.copy()

    # Use numpy vectorized operations
    obstype = np.where(swcm == 5, 247, obstype)  # WVCA/DL
    obstype = np.where(swcm == 3, 246, obstype)  # WVCT
    obstype = np.where(swcm == 2, 251, obstype)  # VIS
    obstype = np.where(swcm == 1, 245, obstype)  # IRLW

    condition = np.logical_and(swcm == 1, chanfreq >= 5e13)  # IRSW
    obstype = np.where(condition, 240, obstype)

    if not np.any(np.isin(obstype, [247, 246, 251, 245, 240])):
        raise ValueError("Error: Unassigned ObsType found ... ")

    return obstype

def _make_obs(comm, input_path, mapping_path):

    # Get container from mapping file first
    logging(comm, 'INFO', 'Get container from bufr')
    container = bufr.Parser(input_path, mapping_path).parse(comm)

    logging(comm, 'DEBUG', f'container list (original): {container.list()}')
    logging(comm, 'DEBUG', f'all_sub_categories =  {container.all_sub_categories()}')
    logging(comm, 'DEBUG', f'category map =  {container.get_category_map()}')

    # Add new/derived data into container
    for cat in container.all_sub_categories():  

        logging(comm, 'DEBUG', f'category = {cat}')

        satid = container.get('variables/satelliteId', cat)
        if satid.size == 0:
            logging(comm, 'WARNING', f'category {cat[0]} does not exist in input file')
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

            logging(comm, 'DEBUG', f'swcm min/max = {swcm.min()} {swcm.max()}')
            logging(comm, 'DEBUG', f'chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')

            obstype = _get_obs_type(swcm, chanfreq)

            logging(comm, 'DEBUG', f'obstype = {obstype}')
            logging(comm, 'DEBUG', f'obstype min/max =  {obstype.min()} {obstype.max()}')

            paths = container.get_paths('variables/windComputationMethod', cat)
            container.add('variables/obstype_uwind', obstype, paths, cat)
            container.add('variables/obstype_vwind', obstype, paths, cat)

            # Add new variables: ObsValue/windEastward & ObsValue/windNorthward 
            wdir = container.get('variables/windDirection', cat)
            wspd = container.get('variables/windSpeed', cat)

            logging(comm, 'DEBUG', f'wdir min/max = {wdir.min()} {wdir.max()}')
            logging(comm, 'DEBUG', f'wspd min/max = {wspd.min()} {wspd.max()}')

            uob, vob = compute_wind_components(wdir, wspd)

            logging(comm, 'DEBUG', f'uob min/max = {uob.min()} {uob.max()}')
            logging(comm, 'DEBUG', f'vob min/max = {vob.min()} {vob.max()}')

            paths = container.get_paths('variables/windSpeed', cat)
            container.add('variables/windEastward', uob, paths, cat)
            container.add('variables/windNorthward', vob, paths, cat)

    # Check
    logging(comm, 'DEBUG', f'container list (updated): {container.list()}')
    logging(comm, 'DEBUG', f'all_sub_categories {container.all_sub_categories()}')

    return container

def create_obs_group(input_path, mapping_path, category, env):

    comm = bufr.mpi.Comm(env["comm_name"])

    description = _make_description(mapping_path, update=True)

    # Check the cache for the data and return it if it exists
    logging(comm, 'DEBUG', f'Check if bufr.DataCache exists? {bufr.DataCache.has(input_path, mapping_path)}')
    if bufr.DataCache.has(input_path, mapping_path):
        container = bufr.DataCache.get(input_path, mapping_path)
        logging(comm, 'INFO', f'Encode {category} from cache')
        data = iodaEncoder(description).encode(container)[(category,)]
        logging(comm, 'INFO', f'Mark {category} as finished in the cache')
        bufr.DataCache.mark_finished(input_path, mapping_path, [category])
        logging(comm, 'INFO', f'Return the encoded data for {category}')
        return data

    container = _make_obs(comm, input_path, mapping_path)

    # Gather data from all tasks into all tasks. Each task will have the complete record 
    logging(comm, 'INFO', f'Gather data from all tasks into all tasks')
    container.all_gather(comm)

    logging(comm, 'INFO', f'Add container to cache')
    # Add the container to the cache
    bufr.DataCache.add(input_path, mapping_path, container.all_sub_categories(), container)

    # Encode the data
    logging(comm, 'INFO', f'Encode {category}')
    data = iodaEncoder(description).encode(container)[(category,)]

    logging(comm, 'INFO', f'Mark {category} as finished in the cache')
    # Mark the data as finished in the cache
    bufr.DataCache.mark_finished(input_path, mapping_path, [category])

    logging(comm, 'INFO', f'Return the encoded data for {category}')
    return data

def create_obs_file(input_path, mapping_path, output_path):

    comm = bufr.mpi.Comm("world")
    container = _make_obs(comm, input_path, mapping_path)
    container.gather(comm)

    description = _make_description(mapping_path, update=True)

    # Encode the data
    if comm.rank() == 0:
        netcdfEncoder(description).encode(container, output_path) 

    logging(comm, 'INFO', f'Return the encoded data')


if __name__ == '__main__':

    start_time = time.time()

    bufr.mpi.App(sys.argv)
    comm = bufr.mpi.Comm("world")

    # Required input arguments as positional arguments
    parser = argparse.ArgumentParser(description="Convert BUFR to NetCDF using a mapping file.")
    parser.add_argument('input', type=str, help='Input BUFR file')
    parser.add_argument('mapping', type=str, help='BUFR2IODA Mapping File')
    parser.add_argument('output', type=str, help='Output NetCDF file')

    args = parser.parse_args()
    infile = args.input
    mapping = args.mapping
    output = args.output

    create_obs_file(infile, mapping, output)

    end_time = time.time()
    running_time = end_time - start_time
    logging(comm, 'INFO', f'Total running time: {running_time}')
