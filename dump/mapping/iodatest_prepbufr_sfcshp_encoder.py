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


def Compute_sequenceNumber(typ, t29):
    sequenceNumber = np.zeros(typ.shape, dtype=np.int32)
    for i in range(len(typ)):
        if (typ[i] == 180 or typ[i] == 280):
            if (t29[i] > 555 and t29[i] < 565):
                sequenceNumber[i] = 0
            else:
                sequenceNumber[i] = 1

    return sequenceNumber


def create_obs_group(cycle_time,input_mapping,input_path):
    CYCLE_TIME = cycle_time
    YAML_PATH = input_mapping #"./iodatest_prepbufr_sfcshp_mapping.yaml"
    INPUT_PATH = input_path
    print(" CYCLE_TIME: ", CYCLE_TIME)

    container = bufr.Parser(INPUT_PATH, YAML_PATH).parse()

    print(" Do DateTime calculation")
    otmct = container.get('variables/obsTimeMinusCycleTime') #dhr
    otmct_paths = container.get_paths('variables/obsTimeMinusCycleTime')
    otmct2 = np.array(otmct)
    cycleTimeSinceEpoch = np.int64(calendar.timegm(time.strptime(str(int(CYCLE_TIME)), '%Y%m%d%H')))
    dateTime = Compute_dateTime(cycleTimeSinceEpoch, otmct2)

    print(" Do sequenceNumber (Obs SubType) calculation")
    typ = container.get('variables/observationType')
    typ_paths = container.get_paths('variables/observationType')
    t29 = container.get('variables/obssubtype')
    t29_paths = container.get_paths('variables/obssubtype')
    seqNum = Compute_sequenceNumber(typ, t29)

    print(" Do tsen and tv calculation")
    tpc = container.get('variables/temperatureEventCode')
    tob = container.get('variables/airTemperatureObsValue')
    tsen = np.full(tob.shape[0], tob.fill_value)
    tsen = np.where(((tpc >= 1) & (tpc < 8)), tob, tsen)
    tvo = np.full(tob.shape[0], tob.fill_value)
    tvo = np.where((tpc == 8), tob, tvo)

    print(" Do tsen and tv QM calculations")
    tobqm = container.get('variables/airTemperatureQualityMarker')
    tsenqm = np.full(tobqm.shape[0], tobqm.fill_value)
    tsenqm = np.where(((tpc >= 1) & (tpc < 8)), tobqm, tsenqm)
    tvoqm = np.full(tobqm.shape[0], tobqm.fill_value)
    tvoqm = np.where((tpc == 8), tobqm, tvoqm)

    print(" Do tsen and tv ObsError calculations")
    toboe = container.get('variables/airTemperatureObsError')
    tsenoe = np.full(toboe.shape[0], toboe.fill_value)
    tsenoe = np.where(((tpc >= 1) & (tpc < 8)), toboe, tsenoe)
    tvooe = np.full(toboe.shape[0], toboe.fill_value)
    tvooe = np.where((tpc == 8), toboe, tvooe)

    print(" Add variables to container.")
    container.add('variables/dateTime', dateTime, otmct_paths)
    container.add('variables/sensibleTemperatureObsValue', tsen, otmct_paths)
    container.add('variables/virtualTemperatureObsValue', tvo, otmct_paths)
    container.add('variables/sensibleTemperatureQualityMarker', tsenqm, otmct_paths)
    container.add('variables/virtualTemperatureQualityMarker', tvoqm, otmct_paths)
    container.add('variables/sensibleTemperatureObsError', tsenoe, otmct_paths)
    container.add('variables/virtualTemperatureObsError', tvooe, otmct_paths)
    container.add('variables/sequenceNumber', seqNum, typ_paths)

    description = bufr.encoders.Description(YAML_PATH)

    print(" Add container and descriptions to dataset.")
    dataset = next(iter(Encoder(description).encode(container).values()))
    return dataset

