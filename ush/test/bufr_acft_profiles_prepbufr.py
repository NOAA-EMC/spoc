import os
import sys
sys.path.append('/work2/noaa/da/nesposito/backend_20240701/bufr_query/build/lib/python3.10')
sys.path.append('/work2/noaa/da/nesposito/backend_20240701/ioda-bundle/build/lib/python3.10')
sys.path.append('/work2/noaa/da/nesposito/backend_20240701/ioda-bundle/build/lib/python3.10/pyioda')
sys.path.append('/work2/noaa/da/nesposito/backend_20240701/ioda-bundle/build/lib/python3.10/pyiodaconv')
import bufr
from pyioda.ioda.Engines.Bufr import Encoder
import copy
import numpy as np
import numpy.ma as ma
import math
import calendar
import time
from datetime import datetime


def Compute_dateTime(cycleTimeSinceEpoch, dhr):

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


def Compute_typ_other(typ, var):

    typ_var = copy.deepcopy(typ)
    typ_var[(typ_var > 300) & (typ_var < 400)] -= 200
    typ_var[(typ_var > 400) & (typ_var < 500)] -= 300
    typ_var[(typ_var > 500) & (typ_var < 600)] -= 400

    for i in range(len(typ_var)):
        if ma.is_masked(var[i]):
            typ_var[i] = typ_var.fill_value

    return typ_var


def Compute_typ_uv(typ, var):

    typ_var = copy.deepcopy(typ)
    typ_var[(typ_var > 300) & (typ_var < 400)] -= 100
    typ_var[(typ_var > 400) & (typ_var < 500)] -= 200
    typ_var[(typ_var > 500) & (typ_var < 600)] -= 300

    for i in range(len(typ_var)):
        if ma.is_masked(var[i]):
            typ_var[i] = typ_var.fill_value

    return typ_var


def Compute_ialr_if_masked(typ, ialr):

    ialr_bc = copy.deepcopy(ialr)
    for i in range(len(ialr_bc)):
        if ma.is_masked(ialr_bc[i]) and (typ[i] >= 330) and (typ[i] < 340):
            ialr_bc[i] = float(0)

    return ialr_bc


def create_obs_group(cycle_time,input_mapping,input_path):
    CYCLE_TIME = round(cycle_time)
    YAML_PATH = input_mapping #"./iodatest_prepbufr_acft_profiles_mapping.yaml"
    INPUT_PATH = input_path

    container = bufr.Parser(INPUT_PATH, YAML_PATH).parse()

    print(" Do DateTime calculation")
    otmct = container.get('variables/obsTimeMinusCycleTime')
    otmct_paths = container.get_paths('variables/obsTimeMinusCycleTime')
    otmct2 = np.array(otmct)
    cycleTimeSinceEpoch = np.int64(calendar.timegm(time.strptime(str(CYCLE_TIME), '%Y%m%d%H')))
    dateTime = Compute_dateTime(cycleTimeSinceEpoch, otmct2)

    container.add('variables/dateTime', dateTime, otmct_paths)

    print(" Do ObsType calculation")
    ot = container.get('variables/observationType')
    ot_paths = container.get_paths('variables/observationType')

    airTemperature = container.get('variables/airTemperatureObsValue')
    #virtualTemperature = container.get('variables/virtualTemperatureObsValue')
    specificHumidity = container.get('variables/specificHumidityObsValue')
    wind = container.get('variables/windNorthwardObsValue')

    ot_airTemperature = Compute_typ_other(ot, airTemperature) 
    #ot_virtualTemperature = Compute_typ_other(ot, virtualTemperature)
    ot_specificHumidity = Compute_typ_other(ot, specificHumidity)
    ot_wind = Compute_typ_uv(ot, wind)

    print(" Change IALR to 0.0 if masked for bias correction.")
    ialr = container.get('variables/instantaneousAltitudeRate0')
    ialr_paths = container.get_paths('variables/instantaneousAltitudeRate0')
    ialr2 = ma.array(ialr)

    ialr_bc = Compute_ialr_if_masked(ot, ialr2)

    print("Make an array of 0s for MetaData/sequenceNumber")
    sequenceNum = np.zeros(ot.shape, dtype=np.int32)

    print(" Add new variables to container")
    container.add('variables/airTemperatureObservationType', ot_airTemperature, ot_paths)
    #container.add('variables/virtualTemperatureObservationType', ot_virtualTemperature, ot_paths)
    container.add('variables/specificHumidityObservationType', ot_specificHumidity, ot_paths)
    container.add('variables/windObservationType', ot_wind, ot_paths)
    container.add('variables/instantaneousAltitudeRate', ialr_bc, ot_paths) #ialr_paths)
    container.add('variables/sequenceNumber', sequenceNum, ot_paths)

    description = bufr.encoders.Description(YAML_PATH)

    print(" Add container and descriptions to dataset.")
    dataset = next(iter(Encoder(description).encode(container).values()))
    return dataset

