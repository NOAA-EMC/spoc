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

