#!/usr/bin/env python3
import os
import sys
import argparse
import bufr
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
logger = Logger('bufr_aircft_profiles_prepbufr.py', level=log_level, colored_log=False)

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

def _compute_datetime(cycleTimeSinceEpoch, dhr):
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

    return dateTime


def _compute_typ_other(typ, var):
    """
    Compute datatype if the variable is not wind.
    Parameters:
        typ: datatype
        var: obsValue variable
    Returns:
        Masked array of the new datatype
    """

    typ_var = copy.deepcopy(typ)
    typ_var[(typ_var > 300) & (typ_var < 400)] -= 200
    typ_var[(typ_var > 400) & (typ_var < 500)] -= 300
    typ_var[(typ_var > 500) & (typ_var < 600)] -= 400

    for i in range(len(typ_var)):
        if ma.is_masked(var[i]):
            typ_var[i] = typ_var.fill_value

    return typ_var


def _compute_typ_uv(typ, var):
    """
    Compute datatype if the variable is wind.
    Parameters:
        typ: datatype
        var: obsValue variable
    Returns:
        Masked array of the new datatype
    """

    typ_var = copy.deepcopy(typ)
    typ_var[(typ_var > 300) & (typ_var < 400)] -= 100
    typ_var[(typ_var > 400) & (typ_var < 500)] -= 200
    typ_var[(typ_var > 500) & (typ_var < 600)] -= 300

    for i in range(len(typ_var)):
        if ma.is_masked(var[i]):
            typ_var[i] = typ_var.fill_value

    return typ_var


def _compute_ialr_if_masked(typ, ialr):
    """
    Compute instantaneousAltitudeRate (IALR) if it is masked.
    Parameters:
        typ: datatype
        ialr: instantaneousAltitudeRate
    Returns:
        Masked array of the updated instantaneousAltitudeRate
    """

    ialr_bc = copy.deepcopy(ialr)
    for i in range(len(ialr_bc)):
        if ma.is_masked(ialr_bc[i]) and (typ[i] >= 330) and (typ[i] < 340):
            ialr_bc[i] = float(0)

    return ialr_bc


def _make_description(mapping_path, cycle_time, update=False):
    description = bufr.encoders.Description(mapping_path)

    ReferenceTime = np.int64(calendar.timegm(time.strptime(str(int(cycle_time)), '%Y%m%d%H')))

    if update:
        # Define the variables to be added in a list of dictionaries
        variables = [
            {
                'name': 'MetaData/sequenceNumber',
                'source': 'variables/sequenceNumber',
                'units': '1',
                'longName': 'Sequence Number (Obs Subtype)',
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

        #description.add_global(name='datetimeReference', value=str(ReferenceTime))

    return description


def _make_obs(comm, input_path, mapping_path, cycle_time):
    """
    Create the ioda acft_profiles prepbufr observations:
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
    cycle_time: str
            The cycle in YYYYMMDDHH format
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
    logging(comm, 'DEBUG', f'lon update max/min: ${lon.max()}, ${lon.min()}')

    logging(comm, 'DEBUG', f'Do DateTime calculation')
    otmct = container.get('variables/obsTimeMinusCycleTime')
    otmct_paths = container.get_paths('variables/obsTimeMinusCycleTime')
    otmct2 = np.array(otmct)
    cycleTimeSinceEpoch = np.int64(calendar.timegm(time.strptime(str(int(cycle_time)), '%Y%m%d%H')))
    dateTime = _compute_datetime(cycleTimeSinceEpoch, otmct2)
    min_dateTime_ge_zero = min(x for x in dateTime if x > -1)
    logging(comm, 'DEBUG', f'dateTime min/max = {min_dateTime_ge_zero} {dateTime.max()}')

    logging(comm, 'DEBUG', f'Make an array of 0s for MetaData/sequenceNumber')
    sequenceNum = np.zeros(lon.shape, dtype=np.int32)
    logging(comm, 'DEBUG', f'sequenceNum min/max =  {sequenceNum.min()} {sequenceNum.max()}')

    logging(comm, 'DEBUG', f'Compute Obstypes')
    t_ot = container.get('variables/airTemperatureObservationType')
    #tv_ot = container.get('variables/virtualTemperatureObservationType')
    q_ot = container.get('variables/specificHumidityObservationType')
    uv_ot = container.get('variables/windObservationType')
    ot_paths = container.get_paths('variables/airTemperatureObservationType')

    airTemperature = container.get('variables/airTemperatureObsValue')
    #virtualTemperature = container.get('variables/virtualTemperatureObsValue')
    specificHumidity = container.get('variables/specificHumidityObsValue')
    wind = container.get('variables/windNorthwardObsValue')

    ot_airTemperature = _compute_typ_other(t_ot, airTemperature)
    #ot_virtualTemperature = _compute_typ_other(tv_ot, virtualTemperature)
    ot_specificHumidity = _compute_typ_other(q_ot, specificHumidity)
    ot_wind = _compute_typ_uv(uv_ot, wind)

    logging(comm, 'DEBUG', f'Change IALR to 0.0 if masked for bias correction.')
    ialr = container.get('variables/instantaneousAltitudeRate')
    ialr_paths = container.get_paths('variables/instantaneousAltitudeRate')
    ialr2 = ma.array(ialr)

    ialr_bc = _compute_ialr_if_masked(uv_ot, ialr2)

    logging(comm, 'DEBUG', f'Update variables in container')
    container.replace('variables/longitude', lon)
    container.replace('variables/instantaneousAltitudeRate', ialr_bc)#, ot_paths)
    container.replace('variables/airTemperatureObservationType', ot_airTemperature)
    #container.replace('variables/virtualTemperatureObservationType', ot_virtualTemperature)
    container.replace('variables/specificHumidityObservationType', ot_specificHumidity)
    container.replace('variables/windObservationType', ot_wind)

    logging(comm, 'DEBUG', f'Add variables to container')
    container.add('variables/sequenceNumber', sequenceNum, lon_paths)

    # Check
    logging(comm, 'DEBUG', f'container list (updated): {container.list()}')

    return container


def create_obs_group(input_path, mapping_path, cycle_time, env):

    comm = bufr.mpi.Comm(env["comm_name"])

    logging(comm, 'INFO', f'Make description and make obs')
    description = _make_description(mapping_path, cycle_time, update=True)
    container = _make_obs(comm, input_path, mapping_path, cycle_time)

    # Gather data from all tasks into all tasks. Each task will have the complete record
    logging(comm, 'INFO', f'Gather data from all tasks into all tasks')
    container.all_gather(comm)

    logging(comm, 'INFO', f'Encode the data')
    data = next(iter(iodaEncoder(description).encode(container).values()))

    logging(comm, 'INFO', f'Return the encoded data.')
    return data


def create_obs_file(input_path, mapping_path, output_path, cycle_time):

    comm = bufr.mpi.Comm("world")
    container = _make_obs(comm, input_path, mapping_path, cycle_time)
    container.gather(comm)

    description = _make_description(mapping_path, cycle_time, update=True)

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
    parser.add_argument('cycle_time', type=str, help='cycle time in YYYYMMDDHH format')

    args = parser.parse_args()
    infile = args.input
    mapping = args.mapping
    output = args.output
    cycle_time = args.cycle_time

    create_obs_file(infile, mapping, output, cycle_time)

    end_time = time.time()
    running_time = end_time - start_time
    logging(comm, 'INFO', f'Total running time: {running_time}')
