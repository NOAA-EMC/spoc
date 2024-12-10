#!/usr/bin/env python3
import os
import sys
import bufr
from pyioda.ioda.Engines.Bufr import Encoder
import copy
import numpy as np
import numpy.ma as ma
import math
import calendar
import time
from datetime import datetime
from wxflow import Logger

# Initialize Logger
# Get log level from the environment variable, default to 'INFO it not set
log_level = os.getenv('LOG_LEVEL', 'INFO')
logger = Logger('BUFR2IODA_adpsfc_prepbufr.py', level=log_level, colored_log=False)

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
                'name':  'MetaData/dateTime'
                'source': 'variables/dateTime',
                'units': 'seconds since 1970-01-01T00:00:00Z',
                'longName': 'dateTime',
            },
            {
                'name': 'MetaData/sequenceNumber',
                'source': 'variables/sequenceNumber'
                'units': '1',
                'longName': 'Sequence Number (Obs Subtype)',
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


def Compute_dateTime(cycleTimeSinceEpoch, dhr):
    """
    Compute dateTime using the cycleTimeSinceEpoch and Cycle Time
        minus Cycle Time 

    Parameters:
        cycleTimeSinceEpoch: Time of cycle in Epoch Time
        dhr: Observation Time Minus Cycle Time

    Returns:
        Masked array of dateTime values
    """

    int64_fill_value = np.int64(0)

    dateTime = np.zeros(dhr.shape, dtype=np.int64)
    for i in range(len(dateTime)):
        if ma.is_masked(dhr[i]):
            continue
        else:
            dateTime[i] = np.int64(dhr[i]*3600) + cycleTimeSinceEpoch

    dateTime = ma.array(dateTime)
    dateTime = ma.masked_values(dateTime, int64_fill_value)

    return dateTime

def _make_obs(comm, input_path, mapping_path):

    # Get container from mapping file first
    logging(comm, 'INFO', 'Get container from bufr')
    container = bufr.Parser(input_path, mapping_path).parse(comm)

    logging(comm, 'DEBUG', f'container list (original): {container.list()}')
    #logging(comm, 'DEBUG', f'all_sub_categories =  {container.all_sub_categories()}')
    #logging(comm, 'DEBUG', f'category map =  {container.get_category_map()}')

    # Add new/derived data into container
    #for cat in container.all_sub_categories():

    #    logging(comm, 'DEBUG', f'category = {cat}')

    #    satid = container.get('variables/obsTimeMinusCycleTime', cat)
    #    if satid.size == 0:
    #        logging(comm, 'WARNING', f'category {cat[0]} does not exist in input file')
            #paths = container.get_paths('variables/windComputationMethod', cat)
            #obstype = container.get('variables/windComputationMethod', cat)
            #container.add('variables/obstype_uwind', obstype, paths, cat)
            #container.add('variables/obstype_vwind', obstype, paths, cat)

            #paths = container.get_paths('variables/windSpeed', cat)
            #wob = container.get('variables/windSpeed', cat)
            #container.add('variables/windEastward', wob, paths, cat)
            #container.add('variables/windNorthward', wob, paths, cat)

    #    else:
    print(" Do DateTime calculation")
    otmct = container.get('variables/obsTimeMinusCycleTime')
    otmct_paths = container.get_paths('variables/obsTimeMinusCycleTime')
    otmct2 = np.array(otmct)
    cycleTimeSinceEpoch = np.int64(calendar.timegm(time.strptime(str(int(CYCLE_TIME)), '%Y%m%d%H')))
    dateTime = Compute_dateTime(cycleTimeSinceEpoch, otmct2)
    logging(comm, 'DEBUG', f'dateTime min/max = {dateTime.min()} {dateTime.max()}')

    print("Make an array of 0s for MetaData/sequenceNumber")
    sequenceNum = np.zeros(otmct.shape, dtype=np.int32)
    logging(comm, 'DEBUG', f' sequenceNummin/max =  {sequenceNum.min()} {sequenceNum.max()}')

    print(" Add variables to container.")
    container.add('variables/dateTime', dateTime, otmct_paths)
    container.add('variables/sequenceNumber', sequenceNum, otmct_paths)
    #description = bufr.encoders.Description(YAML_PATH)

            #print(" Add container and descriptions to dataset.")
            #dataset = next(iter(Encoder(description).encode(container).values()))


            # Add new variables: ObsType/windEastward & ObsType/windNorthward
#            swcm = container.get('variables/windComputationMethod', cat)
#            chanfreq = container.get('variables/sensorCentralFrequency', cat)
#
#            logging(comm, 'DEBUG', f'swcm min/max = {swcm.min()} {swcm.max()}')
#            logging(comm, 'DEBUG', f'chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')
#
#            obstype = _get_obs_type(swcm, chanfreq)
#
#            logging(comm, 'DEBUG', f'obstype = {obstype}')
#            logging(comm, 'DEBUG', f'obstype min/max =  {obstype.min()} {obstype.max()}')
#
#            paths = container.get_paths('variables/windComputationMethod', cat)
#            container.add('variables/obstype_uwind', obstype, paths, cat)
#            container.add('variables/obstype_vwind', obstype, paths, cat)
#
#            # Add new variables: ObsValue/windEastward & ObsValue/windNorthward
#            wdir = container.get('variables/windDirection', cat)
#            wspd = container.get('variables/windSpeed', cat)
#
#            logging(comm, 'DEBUG', f'wdir min/max = {wdir.min()} {wdir.max()}')
#            logging(comm, 'DEBUG', f'wspd min/max = {wspd.min()} {wspd.max()}')
#
#            uob, vob = compute_wind_components(wdir, wspd)
#
#            logging(comm, 'DEBUG', f'uob min/max = {uob.min()} {uob.max()}')
#            logging(comm, 'DEBUG', f'vob min/max = {vob.min()} {vob.max()}')
#
#            paths = container.get_paths('variables/windSpeed', cat)
#            container.add('variables/windEastward', uob, paths, cat)
#            container.add('variables/windNorthward', vob, paths, cat)

    # Check
    logging(comm, 'DEBUG', f'container list (updated): {container.list()}')
    #logging(comm, 'DEBUG', f'all_sub_categories {container.all_sub_categories()}')

    return container

def create_obs_group(cycle_time,input_mapping,input_path):
    CYCLE_TIME = cycle_time
    YAML_PATH = input_mapping #"./iodatest_prepbufr_adpsfc_mapping.yaml"
    INPUT_PATH = input_path
    print(" CYCLE_TIME: ", CYCLE_TIME)

    container = bufr.Parser(INPUT_PATH, YAML_PATH).parse()

    print(" Do DateTime calculation")
    otmct = container.get('variables/obsTimeMinusCycleTime')
    otmct_paths = container.get_paths('variables/obsTimeMinusCycleTime')
    otmct2 = np.array(otmct)
    cycleTimeSinceEpoch = np.int64(calendar.timegm(time.strptime(str(int(CYCLE_TIME)), '%Y%m%d%H')))
    dateTime = Compute_dateTime(cycleTimeSinceEpoch, otmct2)

    print("Make an array of 0s for MetaData/sequenceNumber")
    sequenceNum = np.zeros(otmct.shape, dtype=np.int32)

    print(" Add variables to container.")
    container.add('variables/dateTime', dateTime, otmct_paths)
    container.add('variables/sequenceNumber', sequenceNum, otmct_paths)
    description = bufr.encoders.Description(YAML_PATH)

    print(" Add container and descriptions to dataset.")
    dataset = next(iter(Encoder(description).encode(container).values()))
    return dataset

#def create_obs_group(input_path, mapping_path, category, env):
def create_obs_group(input_path, mapping_path, env):

    comm = bufr.mpi.Comm(env["comm_name"])

    description = _make_description(mapping_path, update=True)

    # Check the cache for the data and return it if it exists
    logging(comm, 'DEBUG', f'Check if bufr.DataCache exists? {bufr.DataCache.has(input_path, mapping_path)}')
    if bufr.DataCache.has(input_path, mapping_path):
        container = bufr.DataCache.get(input_path, mapping_path)
        #logging(comm, 'INFO', f'Encode {category} from cache')
        logging(comm, 'INFO', f'Encode from cache')
        #data = iodaEncoder(description).encode(container)[(category,)]
        data = iodaEncoder(description).encode(container)
        #logging(comm, 'INFO', f'Mark {category} as finished in the cache')
        logging(comm, 'INFO', f'Mark as finished in the cache')
        #bufr.DataCache.mark_finished(input_path, mapping_path, [category])
        bufr.DataCache.mark_finished(input_path, mapping_path)
        #logging(comm, 'INFO', f'Return the encoded data for {category}')
        logging(comm, 'INFO', f'Return the encoded data.')
        return data

    container = _make_obs(comm, input_path, mapping_path)

    # Gather data from all tasks into all tasks. Each task will have the complete record 
    logging(comm, 'INFO', f'Gather data from all tasks into all tasks')
    container.all_gather(comm)

    logging(comm, 'INFO', f'Add container to cache')
    # Add the container to the cache
    #bufr.DataCache.add(input_path, mapping_path, container.all_sub_categories(), container)
    bufr.DataCache.add(input_path, mapping_path, container)

    # Encode the data
    logging(comm, 'INFO', f'Encode {category}')
    #data = iodaEncoder(description).encode(container)[(category,)]
    data = iodaEncoder(description).encode(container)

    logging(comm, 'INFO', f'Mark {category} as finished in the cache')
    # Mark the data as finished in the cache
    #bufr.DataCache.mark_finished(input_path, mapping_path, [category])
    bufr.DataCache.mark_finished(input_path, mapping_path)

    logging(comm, 'INFO', f'Return the encoded data.')
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

    # Required input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input BUFR', required=True)
    parser.add_argument('-m', '--mapping', type=str, help='BUFR2IODA Mapping File', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output NetCDF', required=True)

    args = parser.parse_args()
    mapping = args.mapping
    infile = args.input
    output = args.output

    create_obs_file(infile, mapping, output)

    end_time = time.time()
    running_time = end_time - start_time
    logging(comm, 'INFO', f'Total running time: {running_time}')
