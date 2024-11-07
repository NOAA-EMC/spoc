import os
import sys
sys.path.append('/work2/noaa/da/pkumar/HERCULES/EMC-bufr-query/bufr-query/build/lib')
sys.path.append('/work2/noaa/da/pkumar/HERCULES/EMC-bufr-query/bufr-query/build/lib/python3.10')
sys.path.append('/work2/noaa/da/pkumar/HERCULES/EMC-bufr-query/bufr-query/build/lib/python3.10/site-packages')

import bufr
from bufr.encoders import netcdf
import copy
import numpy as np
import numpy.ma as ma
import math
import calendar
import time
from datetime import datetime
import argparse


def Mask_typ_for_var(typ, var):

    typ_var = copy.deepcopy(typ)
    for i in range(len(typ_var)):
        if ma.is_masked(var[i]):
            typ_var[i] = typ_var.fill_value

    return typ_var


def create_obs_group(input_path, output_path, mapping_path):

    container = bufr.Parser(input_path, mapping_path).parse()

    print("Get the Variables/Paths")
    lat_path = container.get_paths('variables/latitude')
    typ_path = container.get_paths('variables/observationType')
    typ = container.get('variables/observationType')
    cat = container.get('variables/prepbufrDataLevelCategory')
    tpc = container.get('variables/temperatureEventProgramCode')
    tob = container.get('variables/airTemperature')
    pob = container.get('variables/pressure')
    #tvo = container.get('variables/virtualTemperature')
    uob = container.get('variables/windEastward')
    vob = container.get('variables/windNorthward')
    qob = container.get('variables/specificHumidity')
    pqm = container.get('variables/pressureQualityMarker')
    tqm = container.get('variables/airTemperatureQualityMarker')
    poe = container.get('variables/pressureError')
    toe = container.get('variables/airTemperatureError')

    print("ObsValue- Derived Variable Calculation")
    ps = np.full(pob.shape[0], pob.fill_value)
    ps = np.where(cat == 0, pob, ps)

    tsen = np.full(tob.shape[0], tob.fill_value)
    tsen = np.where(((tpc >= 1) & (tpc < 8)), tob, tsen)

    tvo = np.full(tob.shape[0], tob.fill_value)
    tvo = np.where((tpc == 8), tob, tvo)

    print("QualityMarker- Derived Variable Calculation")
    psqm = np.full(pqm.shape[0], pqm.fill_value)
    psqm = np.where(cat == 0, pqm, psqm)

    tsenqm = np.full(tqm.shape[0], tqm.fill_value)
    tsenqm = np.where(((tpc >= 1) & (tpc < 8)), tqm, tsenqm)

    tvoqm = np.full(tqm.shape[0], tqm.fill_value)
    tvoqm = np.where((tpc == 8), tqm, tvoqm)

    print("ObsError- Derived Variable Calculation")
    psoe = ma.array(np.full(poe.shape[0], poe.fill_value))
    psoe = ma.where(cat == 0, poe, psoe)

    tsenoe = np.full(toe.shape[0], toe.fill_value)
    tsenoe = np.where(((tpc >= 1) & (tpc < 8)), toe, tsenoe)

    tvooe = np.full(toe.shape[0], toe.fill_value)
    tvooe = np.where((tpc == 8), toe, tvooe)

    print("Add Variables to Container")
    container.add('variables/stationPressure', ps, lat_path)
    container.add('variables/stationPressureQualityMarker', psqm, lat_path)
    container.add('variables/stationPressureError', psoe, lat_path)

    container.add('variables/sensibleTemperature', tsen, lat_path)
    container.add('variables/sensibleTemperatureQualityMarker', tsenqm, lat_path)
    container.add('variables/sensibleTemperatureError', tsenoe, lat_path)

    #container.add('variables/virtualTemperature', tvo, lat_path)
    container.add('variables/virtualTemperatureQualityMarker', tvoqm, lat_path)
    container.add('variables/virtualTemperatureError', tvooe, lat_path)

    print("Create ObsType Variables")
    typ_pob = Mask_typ_for_var(typ, pob)
    typ_ps  = Mask_typ_for_var(typ, ps)
    typ_tob = Mask_typ_for_var(typ, tob)
    typ_tvo = Mask_typ_for_var(typ, tvo)
    typ_uob = Mask_typ_for_var(typ, uob)
    typ_vob = Mask_typ_for_var(typ, vob)
    typ_qob = Mask_typ_for_var(typ, qob)

    print("Add ObsType Variables to Container")
    container.add('variables/pressureObsType', typ_pob, typ_path)
    container.add('variables/stationPressureObsType', typ_ps, typ_path)
    container.add('variables/airTemperatureObsType', typ_tob, typ_path)
    container.add('variables/virtualTemperatureObsType', typ_tvo, typ_path)
    container.add('variables/windEastwardObsType', typ_uob, typ_path)
    container.add('variables/windNorthwardObsType', typ_vob, typ_path)
    container.add('variables/specificHumidityObsType', typ_qob, typ_path)

    netcdf.Encoder(mapping_path).encode(container, output_path)


if __name__ == '__main__':

    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input BUFR', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output NetCDF', required=True)
    parser.add_argument('-m', '--mapping', type=str, help='Mapping YAML', required=True)
    args = parser.parse_args()

    infile = args.input
    outfile = args.output
    mappingfile = args.mapping

    create_obs_group(infile, outfile, mappingfile)

    end_time = time.time()
    running_time = end_time - start_time
    print("Total running time: ", running_time)
