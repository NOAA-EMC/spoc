#!/usr/bin/env python3
import os
import sys
import bufr
import argparse
import copy
import numpy as np
import numpy.ma as ma
import math
import calendar
import time
from datetime import datetime
from pyioda.ioda.Engines.Bufr import Encoder as iodaEncoder
from bufr.encoders.netcdf import Encoder as netcdfEncoder
from wxflow import Logger


# Initialize Logger
# Get log level from the environment variable, default to 'INFO it not set
log_level = os.getenv('LOG_LEVEL', 'INFO')
logger = Logger('bufr_sfcshp_prepbufr.py', level=log_level, colored_log=False)

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
                'name': 'MetaData/sequenceNumber',
                'source': 'variables/sequenceNumber',
                'units': '1',
                'longName': 'Sequence Number (Obs Subtype)',
            },
    #        {
    #            'name': 'QualityMarker/airTemperature',
    #            'source': 'variables/sensibleTemperatureQualityMarker',
    #            'units': '1',
    #            'longName': 'Air Temperature Quality Marker',
    #        },
    #        {
    #            'name': 'QualityMarker/virtualTemperature',
    #            'source': 'variables/virtualTemperatureQualityMarker',
    #            'units': '1',
    #            'longName': 'Virtual Temperature Quality Marker',
    #        },
    #        {
    #            'name': 'ObsValue/airTemperature',
    #            'source': 'variables/sensibleTemperatureObsValue',
    #            'longName': 'Air Temperature',
    #            'units': 'K',
    #        },
    #        {
    #            'name': 'ObsValue/virtualTemperature',
    #            'source': 'variables/virtualTemperatureObsValue',
    #            'longName': 'Virtual Temperature',
    #            'units': 'K',
    #        },
    #        {
    #            'name': 'ObsError/airTemperature',
    #            'source': 'variables/sensibleTemperatureObsError',
    #            'longName': 'Temperature Error',
    #            'units': 'K',
    #        },
    #        {
    #            'name': 'ObsError/virtualTemperature',
    #            'source': 'variables/virtualTemperatureObsError',
    #            'longName': 'Temperature Error',
    #            'units': 'K',
    #        }
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


def Compute_sequenceNumber(typ, t29):
    sequenceNumber = np.zeros(typ.shape, dtype=np.int32)
    for i in range(len(typ)):
        if (typ[i] == 180 or typ[i] == 280):
            if (t29[i] > 555 and t29[i] < 565):
                sequenceNumber[i] = 0
            else:
                sequenceNumber[i] = 1

    return sequenceNumber


def _make_obs(comm, input_path, mapping_path):
    """
    Create the ioda sfcshp prepbufr observations:
    - reads values 
    - adds sequenceNum 

    Parameters
    ----------
    comm: object
            The communicator object (e.g., MPI) 
    input_path: str
            The input bufr file
    mapping_path: str
            The input bufr2ioda mapping file
    """

    # Get container from mapping file first
    logging(comm, 'INFO', 'Get container from bufr')
    container = bufr.Parser(input_path, mapping_path).parse(comm)

    logging(comm, 'DEBUG', f'container list (original): {container.list()}')
    logging(comm, 'DEBUG', f'Change longitude range from [0,360] to [-180,180]')
    lon = container.get('variables/longitude')
    lon_paths = container.get_paths('variables/longitude')
    lon[lon>180] -= 360
    lon = ma.round(lon, decimals=2)

    logging(comm, 'DEBUG', f'Do sequenceNumber (Obs SubType) calculation')
    typ = container.get('variables/observationType')
    logging(comm, 'DEBUG', f'2Do sequenceNumber (Obs SubType) calculation')
    typ_paths = container.get_paths('variables/observationType')
    logging(comm, 'DEBUG', f'3Do sequenceNumber (Obs SubType) calculation')
    t29 = container.get('variables/obssubtype')
    logging(comm, 'DEBUG', f'4Do sequenceNumber (Obs SubType) calculation')
    t29_paths = container.get_paths('variables/obssubtype')
    logging(comm, 'DEBUG', f'5Do sequenceNumber (Obs SubType) calculation')
    seqNum = Compute_sequenceNumber(typ, t29)
    logging(comm, 'DEBUG', f' sequenceNum min/max =  {seqNum.min()} {seqNum.max()}')

    logging(comm, 'DEBUG', f'Do tsen and tv calculation')
    tpc = container.get('variables/temperatureEventCode')
    tob = container.get('variables/airTemperatureObsValue')
    tob_paths = container.get_paths('variables/airTemperatureObsValue')
    tsen = np.full(tob.shape[0], tob.fill_value)
    tsen = np.where(((tpc >= 1) & (tpc < 8)), tob, tsen)
    tvo = np.full(tob.shape[0], tob.fill_value)
    tvo = np.where((tpc == 8), tob, tvo)

    logging(comm, 'DEBUG', f'Do tsen and tv QM calculations')
    tobqm = container.get('variables/airTemperatureQualityMarker')
    tsenqm = np.full(tobqm.shape[0], tobqm.fill_value)
    tsenqm = np.where(((tpc >= 1) & (tpc < 8)), tobqm, tsenqm)
    tvoqm = np.full(tobqm.shape[0], tobqm.fill_value)
    tvoqm = np.where((tpc == 8), tobqm, tvoqm)

    logging(comm, 'DEBUG', f'Do tsen and tv ObsError calculations')
    toboe = container.get('variables/airTemperatureObsError')
    tsenoe = np.full(toboe.shape[0], toboe.fill_value)
    tsenoe = np.where(((tpc >= 1) & (tpc < 8)), toboe, tsenoe)
    tvooe = np.full(toboe.shape[0], toboe.fill_value)
    tvooe = np.where((tpc == 8), toboe, tvooe)

    logging(comm, 'DEBUG', f'Update variables in container')
    container.replace('variables/longitude', lon)
    container.replace('variables/airTemperatureObsValue', tsen)
    container.replace('variables/airTemperatureQualityMarker', tsenqm)
    container.replace('variables/airTemperatureObsError', tsenoe)
    container.replace('variables/virtualTemperatureObsValue', tvo)
    container.replace('variables/virtualTemperatureQualityMarker', tvoqm)
    container.replace('variables/virtualTemperatureObsError', tvooe)

    logging(comm, 'DEBUG', f'Add variables to container')
    container.add('variables/sequenceNumber', seqNum, typ_paths)
    #container.add('variables/sensibleTemperatureObsValue', tsen, tob_paths)
    #container.add('variables/virtualTemperatureObsValue', tvo, tob_paths)
    #container.add('variables/sensibleTemperatureQualityMarker', tsenqm, tob_paths)
    #container.add('variables/virtualTemperatureQualityMarker', tvoqm, tob_paths)
    #container.add('variables/sensibleTemperatureObsError', tsenoe, tob_paths)
    #container.add('variables/virtualTemperatureObsError', tvooe, tob_paths)

    # Check
    logging(comm, 'DEBUG', f'container list (updated): {container.list()}')

    return container



#def create_obs_group(cycle_time,input_mapping,input_path):
#    CYCLE_TIME = cycle_time
#    YAML_PATH = input_mapping #"./iodatest_prepbufr_sfcshp_mapping.yaml"
#    INPUT_PATH = input_path
#    print(" CYCLE_TIME: ", CYCLE_TIME)
#
#    container = bufr.Parser(INPUT_PATH, YAML_PATH).parse()
#
##    print(" Do DateTime calculation")
##    lon = container.get('variables/obsTimeMinusCycleTime') #dhr
##    lon_paths = container.get_paths('variables/obsTimeMinusCycleTime')
##    lon2 = np.array(lon)
##    cycleTimeSinceEpoch = np.int64(calendar.timegm(time.strptime(str(int(CYCLE_TIME)), '%Y%m%d%H')))
##    dateTime = Compute_dateTime(cycleTimeSinceEpoch, lon2)
#
#    print(" Do sequenceNumber (Obs SubType) calculation")
#    typ = container.get('variables/observationType')
#    typ_paths = container.get_paths('variables/observationType')
#    t29 = container.get('variables/obssubtype')
#    t29_paths = container.get_paths('variables/obssubtype')
#    seqNum = Compute_sequenceNumber(typ, t29)
#
#    print(" Do tsen and tv calculation")
#    tpc = container.get('variables/temperatureEventCode')
#    tob = container.get('variables/airTemperatureObsValue')
#    tsen = np.full(tob.shape[0], tob.fill_value)
#    tsen = np.where(((tpc >= 1) & (tpc < 8)), tob, tsen)
#    tvo = np.full(tob.shape[0], tob.fill_value)
#    tvo = np.where((tpc == 8), tob, tvo)
#
#    print(" Do tsen and tv QM calculations")
#    tobqm = container.get('variables/airTemperatureQualityMarker')
#    tsenqm = np.full(tobqm.shape[0], tobqm.fill_value)
#    tsenqm = np.where(((tpc >= 1) & (tpc < 8)), tobqm, tsenqm)
#    tvoqm = np.full(tobqm.shape[0], tobqm.fill_value)
#    tvoqm = np.where((tpc == 8), tobqm, tvoqm)
#
#    print(" Do tsen and tv ObsError calculations")
#    toboe = container.get('variables/airTemperatureObsError')
#    tsenoe = np.full(toboe.shape[0], toboe.fill_value)
#    tsenoe = np.where(((tpc >= 1) & (tpc < 8)), toboe, tsenoe)
#    tvooe = np.full(toboe.shape[0], toboe.fill_value)
#    tvooe = np.where((tpc == 8), toboe, tvooe)
#
#    print(" Add variables to container.")
#    container.add('variables/dateTime', dateTime, lon_paths)
#    container.add('variables/sensibleTemperatureObsValue', tsen, lon_paths)
#    container.add('variables/virtualTemperatureObsValue', tvo, lon_paths)
#    container.add('variables/sensibleTemperatureQualityMarker', tsenqm, lon_paths)
#    container.add('variables/virtualTemperatureQualityMarker', tvoqm, lon_paths)
#    container.add('variables/sensibleTemperatureObsError', tsenoe, lon_paths)
#    container.add('variables/virtualTemperatureObsError', tvooe, lon_paths)
#    container.add('variables/sequenceNumber', seqNum, typ_paths)
#
#    description = bufr.encoders.Description(YAML_PATH)
#
#    print(" Add container and descriptions to dataset.")
#    dataset = next(iter(Encoder(description).encode(container).values()))
#    return dataset

def create_obs_group(input_path, mapping_path, env):

    comm = bufr.mpi.Comm(env["comm_name"])

    logging(comm, 'INFO', f'Make description and make obs')
    description = _make_description(mapping_path, update=True)
    container = _make_obs(comm, input_path, mapping_path)

    # Gather data from all tasks into all tasks. Each task will have the complete record
    logging(comm, 'INFO', f'Gather data from all tasks into all tasks')
    container.all_gather(comm)

    logging(comm, 'INFO', f'Encode the data')
    data = next(iter(iodaEncoder(description).encode(container).values()))

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
