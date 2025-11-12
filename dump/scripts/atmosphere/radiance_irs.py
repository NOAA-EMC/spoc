#!/usr/bin/env python3

import os
import bufr
import numpy as np
from netCDF4 import Dataset
from bufr.encoders import netcdf

dir   = "/scratch3/NCEPDEV/global/Jack.Woollen/IRSPP/IRSPP/IRSPPv1.3_test_cases/input"
out   = "out.nc"
yaml  = 'yaml.h'
first = 'true'
iter  = 1
imax = np.empty(1024,dtype=int)
jmax = np.empty(1024,dtype=int)

for filename in os.listdir(dir):
   print(filename)
   irs  = Dataset(dir+"/"+filename)
   plat = irs['state/platform']
   cele = irs['state/celestial']
   inst = irs['state/instrument']
   loca = irs['data']
   mwmd = irs['data/mwir']
   lwmd = irs['data/lwir']
   mwqa = irs['data/mwir/quality_band']
   lwqa = irs['data/lwir/quality_band']
   mwva = irs['data/mwir/compressed']
   lwva = irs['data/lwir/compressed']

   # select the hottest spot in each 5x5 box within each dwell
   n = 0
   for a in range(0,155,5):
      for b in range(0,155,5):
         pcmax = 0
         for i in range(1,3):
            for j in range(1,3):
               pc1=lwva.variables['global_pc_scores'][i+a][j+b][1]
               pc1=abs(pc1)
               pcmax = max(pc1,pcmax)
               if pcmax == pc1:
                  imax[n]=i+a
                  jmax[n]=j+b
         n=n+1
         #print(n,imax[n],jmax[n])

   # create numpy for this dwell group
   n=n-1
   atime=np.empty(n)
   adwell_number=np.empty(n)
   astroke_direction=np.empty(n)
   alatitude=np.empty(n)
   alongitude=np.empty(n)
   asatellite_azimuth_angle=np.empty(n)
   asatellite_zenith_angle=np.empty(n)
   asolar_azimuth_angle=np.empty(n)
   asolar_zenith_angle=np.empty(n)
   acloud_signal=np.empty(n)
   acloud_fraction=np.empty(n)
   amwir_global_pc_scores=np.empty((n,150))
   amwir_global_pcr_scores=np.empty(n)
   amwir_global_pcrs_quality=np.empty(n)
   amwir_spatial_sample_quality=np.empty(n)
   amwir_residual_energy=np.empty(n)
   alwir_global_pc_scores=np.empty((n,150))
   alwir_global_pcr_scores=np.empty(n)
   alwir_global_pcrs_quality=np.empty(n)
   alwir_spatial_sample_quality=np.empty(n)
   alwir_residual_energy=np.empty(n)

   # save the soundings selected from this dwell
   for m in range(1,n):
      i = imax[m]
      j = jmax[m]
      atime[m]=loca.variables['time'][:]
      adwell_number[m]=loca.variables['dwell_number'][:]
      astroke_direction[m]=loca.variables['stroke_direction'][:]
      alatitude[m]=loca.variables['latitude'][i][j]
      alongitude[m]=loca.variables['longitude'][i][j]
      asatellite_azimuth_angle[m]=loca.variables['satellite_azimuth_angle'][i][j]
      asatellite_zenith_angle[m]=loca.variables['satellite_zenith_angle'][i][j]
      asolar_azimuth_angle[m]=loca.variables['solar_azimuth_angle'][i][j]
      asolar_zenith_angle[m]=loca.variables['solar_zenith_angle'][i][j]
      acloud_signal[m]=loca.variables['cloud_signal'][i][j]
      acloud_fraction[m]=loca.variables['cloud_fraction'][i][j]
      amwir_global_pc_scores[m]=mwva.variables['global_pc_scores'][i][j][:]
      amwir_global_pcr_scores[m]=mwva.variables['global_pcr_scores'][i][j]
      amwir_global_pcrs_quality[m]=mwva.variables['global_pcrs_quality'][i][j]
      amwir_spatial_sample_quality[m]=mwva.variables['spatial_sample_quality'][i][j]
      amwir_residual_energy[m]=mwva.variables['residual_energy'][:]
      alwir_global_pc_scores[m]=lwva.variables['global_pc_scores'][i][j][:]
      alwir_global_pcr_scores[m]=lwva.variables['global_pcr_scores'][i][j]
      alwir_global_pcrs_quality[m]=lwva.variables['global_pcrs_quality'][i][j]
      alwir_spatial_sample_quality[m]=lwva.variables['spatial_sample_quality'][i][j]
      alwir_residual_energy[m]=lwva.variables['residual_energy'][:]

   # accumulate this dwell into the dump group
   if first == 'true':
      first = 'nottrue'
      time=atime
      dwell_number=adwell_number
      stroke_direction=astroke_direction
      latitude=alatitude
      longitude=alongitude
      satellite_azimuth_angle=asatellite_azimuth_angle
      satellite_zenith_angle=asatellite_zenith_angle
      solar_azimuth_angle=asolar_azimuth_angle
      solar_zenith_angle=asolar_zenith_angle
      cloud_signal=acloud_signal
      cloud_fraction=acloud_fraction
      mwir_global_pc_scores=amwir_global_pc_scores
      mwir_global_pcr_scores=amwir_global_pcr_scores
      mwir_global_pcrs_quality=amwir_global_pcrs_quality
      mwir_spatial_sample_quality=amwir_spatial_sample_quality
      mwir_residual_energy=amwir_residual_energy
      lwir_global_pc_scores=alwir_global_pc_scores
      lwir_global_pcr_scores=alwir_global_pcr_scores
      lwir_global_pcrs_quality=alwir_global_pcrs_quality
      lwir_spatial_sample_quality=alwir_spatial_sample_quality
      lwir_residual_energy=alwir_residual_energy
   else:
      time=np.concatenate((time,atime))
      dwell_number=np.concatenate((dwell_number,adwell_number))
      stroke_direction=np.concatenate((stroke_direction,astroke_direction))
      latitude=np.concatenate((latitude,alatitude))
      longitude=np.concatenate((longitude,alongitude))
      satellite_azimuth_angle=np.concatenate((satellite_azimuth_angle,asatellite_azimuth_angle))
      satellite_zenith_angle=np.concatenate((satellite_zenith_angle,asatellite_zenith_angle))
      solar_azimuth_angle=np.concatenate((solar_azimuth_angle,asolar_azimuth_angle))
      solar_zenith_angle=np.concatenate((solar_zenith_angle,asolar_zenith_angle))
      cloud_signal=np.concatenate((cloud_signal,acloud_signal))
      cloud_fraction=np.concatenate((cloud_fraction,acloud_fraction))
      mwir_global_pc_scores=np.concatenate((mwir_global_pc_scores,amwir_global_pc_scores))
      mwir_global_pcr_scores=np.concatenate((mwir_global_pcr_scores,amwir_global_pcr_scores))
      mwir_global_pcrs_quality=np.concatenate((mwir_global_pcrs_quality,amwir_global_pcrs_quality))
      mwir_spatial_sample_quality=np.concatenate((mwir_spatial_sample_quality,amwir_spatial_sample_quality))
      mwir_residual_energy=np.concatenate((mwir_residual_energy,amwir_residual_energy))
      lwir_global_pc_scores=np.concatenate((lwir_global_pc_scores,alwir_global_pc_scores))
      lwir_global_pcr_scores=np.concatenate((lwir_global_pcr_scores,alwir_global_pcr_scores))
      lwir_global_pcrs_quality=np.concatenate((lwir_global_pcrs_quality,alwir_global_pcrs_quality))
      lwir_spatial_sample_quality=np.concatenate((lwir_spatial_sample_quality,alwir_spatial_sample_quality))
      lwir_residual_energy=np.concatenate((lwir_residual_energy,alwir_residual_energy))

#  if iter==1:
#     break
#  iter=iter+1

# change dtypes for certain variables
stroke_direction = stroke_direction.astype('i')
mwir_global_pcrs_quality = mwir_global_pcrs_quality.astype('i')
mwir_spatial_sample_quality = mwir_spatial_sample_quality.astype('i')
lwir_global_pcrs_quality = lwir_global_pcrs_quality.astype('i')
lwir_spatial_sample_quality = lwir_spatial_sample_quality.astype('i')

# write into the container and the ioda dump file
container = bufr.DataContainer()
description = bufr.encoders.Description(yaml)
container.add('time', time, ['*'])
container.add('dwell_number', dwell_number, ['*'])
container.add('stroke_direction', stroke_direction, ['*'])
container.add('latitude', latitude, ['*'])
container.add('longitude', longitude, ['*'])
container.add('satellite_azimuth_angle', satellite_azimuth_angle, ['*'])
container.add('satellite_zenith_angle', satellite_zenith_angle, ['*'])
container.add('solar_azimuth_angle', solar_azimuth_angle, ['*'])
container.add('solar_zenith_angle', solar_zenith_angle, ['*'])
container.add('cloud_signal', cloud_signal, ['*'])
container.add('cloud_fraction', cloud_fraction, ['*'])
container.add('mwir_global_pc_scores', mwir_global_pc_scores, ['*','*/PSCORE'])
container.add('mwir_global_pcr_scores', mwir_global_pcr_scores, ['*'])
container.add('mwir_global_pcrs_quality', mwir_global_pcrs_quality, ['*'])
container.add('mwir_spatial_sample_quality', mwir_spatial_sample_quality, ['*'])
container.add('mwir_residual_energy', mwir_residual_energy, ['*'])
container.add('lwir_global_pc_scores', lwir_global_pc_scores, ['*','*/PSCORE'])
container.add('lwir_global_pcr_scores', lwir_global_pcr_scores, ['*'])
container.add('lwir_global_pcrs_quality', lwir_global_pcrs_quality, ['*'])
container.add('lwir_spatial_sample_quality', lwir_spatial_sample_quality, ['*'])
container.add('lwir_residual_energy', lwir_residual_energy, ['*'])
netcdf.Encoder(description).encode(container,out)

