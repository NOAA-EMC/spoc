#!/usr/bin/env python3 # read mtg/irs netcdf files compute radiances and write to obsforge ioda file

import os
import sys
import bufr
import rddump
import ncepbufr
import numpy as np
from bufr.encoders import netcdf
from datetime import datetime, timezone

# ----------------------------------------------
# define the input and  parameter file locations
# ----------------------------------------------

filename = 'work/aircar.2023010100'
iodout   = 'aircraft.nc'
yaml     = 'aircraft.yaml'

# -------------------------------------------------
# open the aircraft dump file and count the reports
# -------------------------------------------------

nloc=rddump.rddump(filename)
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
md_observationTypeNum = np.empty(nloc,dtype=int)
md_observationSubTypeNum = np.empty(nloc,dtype=int)
md_latitude = np.empty(nloc,dtype=float)
md_longitude = np.empty(nloc,dtype=float)
md_unix_time = np.empty(nloc,dtype=int)
md_pressure = np.empty(nloc,dtype=float)
md_elevation = np.empty(nloc,dtype=float)
md_instantaneousAltitudeRate = np.empty(nloc,dtype=float)
md_aircraftNavigationSystem = np.empty(nloc,dtype=int)
md_aircraftFlightPhase = np.empty(nloc,dtype=int)

ov_airTemperature = np.empty(nloc,dtype=float)
ov_specificHumidity = np.empty(nloc,dtype=float)
ov_windEastward = np.empty(nloc,dtype=float)
ov_windNorthward = np.empty(nloc,dtype=float)

qm_airTemperature = np.empty(nloc,dtype=float)
qm_specificHumidity = np.empty(nloc,dtype=float)
qm_windEastward = np.empty(nloc,dtype=float)
qm_windNorthward = np.empty(nloc,dtype=float)

oe_airTemperature = np.empty(nloc,dtype=float)
oe_specificHumidity = np.empty(nloc,dtype=float)
oe_windEastward = np.empty(nloc,dtype=float)
oe_windNorthward = np.empty(nloc,dtype=float)

# -------------------------------------------------------------
# define the bufr reader common connection and read in the data 
# -------------------------------------------------------------

common_data = rddump.data

for n in range(nloc):

   next=rddump.rddump('readns')
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
container.add('stationIdentification',      md_stationIdentification,['*'])
container.add('aircraftFlightNumber',       md_aircraftFlightNumber,['*'])
container.add('aircraftTailNumber',         md_aircraftTailNumber,['*'])
container.add('observationTypeNum',         md_observationTypeNum,['*'])
container.add('observationSubTypeNum',      md_observationSubTypeNum,['*'])
container.add('latitude',                   md_latitude,['*'])
container.add('longitude',                  md_longitude,['*'])
container.add('unix_time',                  md_unix_time,['*'])
container.add('pressure',                   md_pressure,['*'])
container.add('elevation',                  md_elevation,['*'])
container.add('aircraftFlightPhase',        md_aircraftFlightPhase,['*'])
container.add('instantaneousAltitudeRate',  md_instantaneousAltitudeRate,['*'])
container.add('aircraftNavigationSystem',   md_aircraftNavigationSystem,['*'])
container.add('airTemperature',             ov_airTemperature,['*'])
container.add('specificHumidity',           ov_specificHumidity,['*'])
container.add('windEastward',               ov_windEastward,['*'])
container.add('windNorthward',              ov_windNorthward,['*'])
container.add('airTemperature',             qm_airTemperature,['*'])
container.add('specificHumidity',           qm_specificHumidity,['*'])
container.add('windEastward',               qm_windEastward,['*'])
container.add('windNorthward',              qm_windNorthward,['*'])
container.add('airTemperature',             oe_airTemperature,['*'])
container.add('specificHumidity',           oe_specificHumidity,['*'])
container.add('windEastward',               oe_windEastward,['*'])
container.add('windNorthward',              oe_windNorthward,['*'])
netcdf.Encoder(description).encode(container,iodout)

