#!/usr/bin/env python3 # read mtg/irs netcdf files compute radiances and write to obsforge ioda file

import os
import sys
import bufr
import ncepbufr
import numpy as np
import dumpbufr_aircraft
from bufr.encoders import netcdf
from datetime import datetime, timezone

# ----------------------------------------------
# define the input and  parameter file locations
# ----------------------------------------------

filename = 'work/acdump.2023080100'
yaml     = '/scratch3/NCEPDEV/global/Jack.Woollen/spoc/dump/config/atmosphere/dumpbufr_aircraft.yaml'
iodout   = 'dumpbufr_aircraft.nc'

# -------------------------------------------------
# open the aircraft dump file and count the reports
# -------------------------------------------------

nloc=dumpbufr_aircraft.rddump(filename)
if nloc>=0:
   print()
   print('Sucessfully opened ',filename,' with ',nloc,' subsets')
   print()
else:
   print('Problem opening ',filename)
   exit()

# ----------------------------------------------------
# define the numpy arrays to capture the aircraft data
# ----------------------------------------------------

md_stationIdentification = np.empty(nloc,dtype='U8')
md_aircraftFlightNumber = np.empty(nloc,dtype='U8')
md_aircraftTailNumber = np.empty(nloc,dtype='U8')
md_observationTypeNum = np.empty(nloc,dtype=np.int32)
md_observationSubTypeNum = np.empty(nloc,dtype=np.int32)
md_latitude = np.empty(nloc,dtype=np.float32)
md_longitude = np.empty(nloc,dtype=np.float32)
md_unix_time = np.empty(nloc,dtype=np.int32)
md_pressure = np.empty(nloc,dtype=np.float32)
md_elevation = np.empty(nloc,dtype=np.float32)
md_instantaneousAltitudeRate = np.empty(nloc,dtype=np.float32)
md_aircraftNavigationSystem = np.empty(nloc,dtype=np.int32)
md_aircraftFlightPhase = np.empty(nloc,dtype=np.int32)

ov_airTemperature = np.empty(nloc,dtype=np.float32)
ov_specificHumidity = np.empty(nloc,dtype=np.float32)
ov_windEastward = np.empty(nloc,dtype=np.float32)
ov_windNorthward = np.empty(nloc,dtype=np.float32)

qm_airTemperature = np.empty(nloc,dtype=np.float32)
qm_specificHumidity = np.empty(nloc,dtype=np.float32)
qm_windEastward = np.empty(nloc,dtype=np.float32)
qm_windNorthward = np.empty(nloc,dtype=np.float32)

oe_airTemperature = np.empty(nloc,dtype=np.float32)
oe_specificHumidity = np.empty(nloc,dtype=np.float32)
oe_windEastward = np.empty(nloc,dtype=np.float32)
oe_windNorthward = np.empty(nloc,dtype=np.float32)

# -------------------------------------------------------------
# define the bufr reader common connection and read in the data 
# -------------------------------------------------------------

common_data = dumpbufr_aircraft.data

for n in range(nloc):

   next=dumpbufr_aircraft.rddump('readns')
   if next != 0:
      print('error reading bufr file')
      sys.exit(99)
   
   year = int(common_data.year)
   mnth = int(common_data.mnth)
   days = int(common_data.days)
   hour = int(common_data.hour)
   minu = int(common_data.minu)
   seco = int(common_data.seco)
   unix_time = datetime(year,mnth,days,hour,minu,seco,tzinfo=timezone.utc).timestamp()

   stid = common_data.stid
   acfn = common_data.acfn 
   actn = common_data.actn

# MetaData

   md_stationIdentification[n]     = str(stid,'utf-8')
   md_aircraftFlightNumber[n]      = str(acfn,'utf-8')
   md_aircraftTailNumber[n]        = str(actn,'utf-8')
   md_observationTypeNum[n]        = int(common_data.otyp)
   md_observationSubTypeNum[n]     = int(common_data.styp)
   md_latitude[n]                  = common_data.flat
   md_longitude[n]                 = common_data.flon
   md_unix_time[n]                 = unix_time
   md_pressure[n]                  = common_data.pres
   md_elevation[n]                 = common_data.elev
   md_aircraftFlightPhase[n]       = common_data.poaf
   md_instantaneousAltitudeRate[n] = common_data.ialr
   md_aircraftNavigationSystem[n]  = common_data.acns

# ObsValue

   ov_airTemperature[n]   = common_data.ovat
   ov_specificHumidity[n] = common_data.ovsh
   ov_windEastward[n]     = common_data.ovew
   ov_windNorthward[n]    = common_data.ovnw

# QualityMarker

   qm_airTemperature[n]   = common_data.qmat
   qm_specificHumidity[n] = common_data.qmsh
   qm_windEastward[n]     = common_data.qmew
   qm_windNorthward[n]    = common_data.qmnw

# ObsError  

   oe_airTemperature[n]   = common_data.oeat
   oe_specificHumidity[n] = common_data.oesh
   oe_windEastward[n]     = common_data.oeew
   oe_windNorthward[n]    = common_data.oenw

# -----------------------------------------
# option to print each report one at a time  
# -----------------------------------------

   if 0==1:
      print('-------------------------------------------------')
      print(md_stationIdentification[n])
      print(md_aircraftFlightNumber[n])
      print(md_aircraftTailNumber[n])
      print(md_observationTypeNum[n])             
      print(md_observationSubTypeNum[n])
      print(md_latitude[n])
      print(md_longitude[n])
      print(md_pressure[n])
      print(md_elevation[n])
      print(md_unix_time[n])
      print(md_aircraftFlightPhase[n]) 
      print(md_instantaneousAltitudeRate[n])
      print(md_aircraftNavigationSystem[n])
      print()
      print(ov_airTemperature[n]) 
      print(ov_specificHumidity[n]) 
      print(ov_windEastward[n]) 
      print(ov_windNorthward[n])
      print()
      print(qm_airTemperature[n]) 
      print(qm_specificHumidity[n]) 
      print(qm_windEastward[n]) 
      print(qm_windNorthward[n])
      print()
      print(oe_airTemperature[n]) 
      print(oe_specificHumidity[n]) 
      print(oe_windEastward[n]) 
      print(oe_windNorthward[n])
      input()

# -------------------------------------------
# build the ioda container and write the file
# -------------------------------------------

container = bufr.DataContainer()
description = bufr.encoders.Description(yaml)
container.add('md_stationIdentification',      md_stationIdentification,['*'])
container.add('md_aircraftFlightNumber',       md_aircraftFlightNumber,['*'])
container.add('md_aircraftTailNumber',         md_aircraftTailNumber,['*'])
container.add('md_observationTypeNum',         md_observationTypeNum,['*'])
container.add('md_observationSubTypeNum',      md_observationSubTypeNum,['*'])
container.add('md_latitude',                   md_latitude,['*'])
container.add('md_longitude',                  md_longitude,['*'])
container.add('md_unix_time',                  md_unix_time,['*'])
container.add('md_pressure',                   md_pressure,['*'])
container.add('md_elevation',                  md_elevation,['*'])
container.add('md_aircraftFlightPhase',        md_aircraftFlightPhase,['*'])
container.add('md_instantaneousAltitudeRate',  md_instantaneousAltitudeRate,['*'])
container.add('md_aircraftNavigationSystem',   md_aircraftNavigationSystem,['*'])
container.add('ov_airTemperature',             ov_airTemperature,['*'])
container.add('ov_specificHumidity',           ov_specificHumidity,['*'])
container.add('ov_windEastward',               ov_windEastward,['*'])
container.add('ov_windNorthward',              ov_windNorthward,['*'])
container.add('qm_airTemperature',             qm_airTemperature,['*'])
container.add('qm_specificHumidity',           qm_specificHumidity,['*'])
container.add('qm_windEastward',               qm_windEastward,['*'])
container.add('qm_windNorthward',              qm_windNorthward,['*'])
container.add('oe_airTemperature',             oe_airTemperature,['*'])
container.add('oe_specificHumidity',           oe_specificHumidity,['*'])
container.add('oe_windEastward',               oe_windEastward,['*'])
container.add('oe_windNorthward',              oe_windNorthward,['*'])
netcdf.Encoder(description).encode(container,iodout)

