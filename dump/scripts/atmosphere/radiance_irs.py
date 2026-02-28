#!/usr/bin/env python3 # read mtg/irs netcdf files compute radiances and write to obsforge ioda file

import os
import sys
import bufr
import time  as tm
import numpy as np
from netCDF4 import Dataset
from bufr.encoders import netcdf

#---define the apodising function

def apply_hamming(x,nchan):
     hamming0 = np.float64(0.54); hamming1 = np.float64(0.23)
     y = x.copy(); ncn = nchan-1
     y[0] = (x[0]*hamming0+x[1]*hamming1)/(hamming0+hamming1)
     y[ncn] = (x[ncn]*hamming0+x[ncn-1]*hamming1)/(hamming0+hamming1)
     for icn in range(1,ncn):
          y[icn] = x[icn-1]*hamming1+x[icn]*hamming0+x[icn+1]*hamming1
     return y

#---read two arguments to define task
#---1) the input data directory path
#---2) the output path/filename

if len(sys.argv) < 3:
    print(f"{sys.argv[0]} needs <inpdir> and <iodout> ")
    exit()

inpdir=sys.argv[1]; print(sys.argv[1])
iodout=sys.argv[2]; print(sys.argv[2])

#---setup more filenames and parameters

yamls   = "/scratch3/NCEPDEV/global/Jack.Woollen/spoc/dump/config/atmosphere"
parms   = "/scratch3/NCEPDEV/global/Jack.Woollen/spoc/dump/scripts/atmosphere/parm"
lwchans = os.path.join(parms,"irs_coopman_lw_channels.txt")
mwchans = os.path.join(parms,"irs_coopman_mw_channels.txt")
reconst = os.path.join(parms,"RSP_OPE_BASEEV_MTS1+IRS_20230925000000_V1_out.h5")
yaml    = os.path.join(yamls,"radiance_irs.yaml")

#---read the channel data for lw and mw ir wave numbers

chan_lw = np.loadtxt(lwchans,dtype=int); chan_lw = chan_lw[chan_lw != 0]; chns_lw = chan_lw.shape[0]
chan_mw = np.loadtxt(mwchans,dtype=int); chan_mw = chan_mw[chan_mw != 0]; chns_mw = chan_mw.shape[0]

#---read reconstruction means and operators for lw and mw ir channels

hd5file=Dataset(reconst)
means = hd5file['lwir'].variables['Mean'][:][:]; means_lw = means[chan_lw]
means = hd5file['mwir'].variables['Mean'][:][:]; means_mw = means[chan_mw]
recop = hd5file['lwir'].variables['ReconstructionOperator']; recop_lw = recop[:,chan_lw]; wnum_lw = recop.shape[1]
recop = hd5file['mwir'].variables['ReconstructionOperator']; recop_mw = recop[:,chan_mw]; wnum_mw = recop.shape[1]

#---hamming the eigenvectors

means_lw = apply_hamming(means_lw,chns_lw)
means_mw = apply_hamming(means_mw,chns_mw)
for ipc in range(len(recop_lw)):
     recop_lw[ipc] = apply_hamming(recop_lw[ipc][:],chns_lw)
     recop_mw[ipc] = apply_hamming(recop_mw[ipc][:],chns_mw)

#---setup radiance channel output parameters 

lw0 = 0
lw1 = chns_lw
mw0 = chns_lw
mw1 = mw0+chns_mw
nchan = chns_lw+chns_mw
chans = np.full(nchan,0,dtype=int)
chans[lw0:lw1] = chan_lw[:]
chans[mw0:mw1] = chan_mw[:]+wnum_lw

#---working parameters and arrays

kmax   = 1024 
jmax   = np.zeros(kmax,dtype=int)
imax   = np.zeros(kmax,dtype=int)
fill   = 1.e300

#---check that the files in the list exist and are readable

try:
   file_list = os.listdir(inpdir)
   file_keep = [] 
except OSError as e:
   sys.stderr.write(f"Error reading input directory '{inpdir}': {e}\n")
   sys.exit(1)

netcdf_files = [f for f in file_list if f.lower().endswith(".nc")]
if not netcdf_files:
   sys.stderr.write(f"No NetCDF files found in input directory '{inpdir}'.\n")
   sys.exit(1)

for filename in netcdf_files:
   try:
        irs  = Dataset(os.path.join(inpdir, filename))
        file_keep.append(filename)
   except OSError as e:
        sys.stderr.write(f"Skipping file '{filename}': failed to open as NetCDF: {e}\n")
        continue

nfil=len(file_keep)      ; print('nfiles=',nfil)

#---create empty numpy array for all the dwell groups

nloc=0 # to start with 

time=np.zeros(nloc)
dwell_number=np.zeros(nloc)
stroke_direction=np.zeros(nloc)
latitude=np.zeros(nloc)
longitude=np.zeros(nloc)
satellite_azimuth_angle=np.zeros(nloc)
satellite_zenith_angle=np.zeros(nloc)
solar_azimuth_angle=np.zeros(nloc)
solar_zenith_angle=np.zeros(nloc)
cloud_signal=np.zeros(nloc)
cloud_fraction=np.zeros(nloc)
mwir_global_pc_scores=np.zeros((nloc,150))
mwir_global_pcr_scores=np.zeros(nloc)
mwir_global_pcrs_quality=np.zeros(nloc)
mwir_spatial_sample_quality=np.zeros(nloc)
mwir_residual_energy=np.zeros(nloc)
lwir_global_pc_scores=np.zeros((nloc,150))
lwir_global_pcr_scores=np.zeros(nloc)
lwir_global_pcrs_quality=np.zeros(nloc)
lwir_spatial_sample_quality=np.zeros(nloc)
lwir_residual_energy=np.zeros(nloc)
radiance=np.zeros((nloc,nchan))
chan_num=np.zeros((nloc,nchan))

#---loop through the list of mtg-irs dwell files to be extracted


dwel=0
for filename in file_keep:
   print(filename)
   irs  = Dataset(os.path.join(inpdir,filename))
   loca = irs['data']
   mwva = irs['data/mwir/compressed']
   lwva = irs['data/lwir/compressed']
   dwel = dwel+1

#---select the hottest spot in each 5x5 box within the dwell

   atime=0; atime=tm.time()

   source = lwva.variables["global_pc_scores"][:][:][:]

   n = 0
   for a in range(0,160,5):
      for b in range(0,160,5):
         pcmax = -99.e99
         for i in range(1,4):
            for j in range(1,4):
               pc1=source[i+a][j+b][0]
               if abs(pc1) < fill: 
                  pcmax = max(pc1,pcmax)
                  if pcmax == pc1:
                     imax[n]=i+a
                     jmax[n]=j+b
         n=n+1

#---save the soundings selected from this dwell

   btime=0; btime=tm.time()

   scalar    = np.zeros(kmax)
   scalar[:] = loca.variables["time"                       ]   [:]; time=np.concatenate((time,scalar))
   scalar[:] = loca.variables["dwell_number"               ]   [:]; dwell_number=np.concatenate((dwell_number,scalar))
   scalar[:] = loca.variables["stroke_direction"           ]   [:]; stroke_direction=np.concatenate((stroke_direction,scalar))
   scalar[:] = mwva.variables["residual_energy"            ]   [:]; mwir_residual_energy=np.concatenate((lwir_residual_energy,scalar))
   scalar[:] = lwva.variables["residual_energy"            ]   [:]; lwir_residual_energy=np.concatenate((lwir_residual_energy,scalar))

   source = loca.variables["latitude"                      ][:][:]; latitude=np.concatenate((latitude,source[imax,jmax]))
   source = loca.variables["longitude"                     ][:][:]; longitude=np.concatenate((longitude,source[imax,jmax]))
   source = loca.variables["satellite_azimuth_angle"       ][:][:]; satellite_azimuth_angle=np.concatenate((satellite_azimuth_angle,source[imax,jmax]))
   source = loca.variables["satellite_zenith_angle"        ][:][:]; satellite_zenith_angle=np.concatenate((satellite_zenith_angle,source[imax,jmax]))
   source = loca.variables["solar_azimuth_angle"           ][:][:]; solar_azimuth_angle=np.concatenate((solar_azimuth_angle,source[imax,jmax]))
   source = loca.variables["solar_zenith_angle"            ][:][:]; solar_zenith_angle=np.concatenate((solar_zenith_angle,source[imax,jmax]))
   source = loca.variables["cloud_signal"                  ][:][:]; cloud_signal=np.concatenate((cloud_signal,source[imax,jmax]))
   source = loca.variables["cloud_fraction"                ][:][:]; cloud_fraction=np.concatenate((cloud_fraction,source[imax,jmax]))
   source = mwva.variables["global_pcr_scores"             ][:][:]; mwir_global_pcr_scores=np.concatenate((mwir_global_pcr_scores,source[imax,jmax]))
   source = mwva.variables["global_pcrs_quality"           ][:][:]; mwir_global_pcrs_quality=np.concatenate((mwir_global_pcrs_quality,source[imax,jmax]))
   source = mwva.variables["spatial_sample_quality"        ][:][:]; mwir_spatial_sample_quality=np.concatenate((mwir_spatial_sample_quality,source[imax,jmax]))
   source = lwva.variables["global_pcr_scores"             ][:][:]; lwir_global_pcr_scores=np.concatenate((lwir_global_pcr_scores,source[imax,jmax]))
   source = lwva.variables["global_pcrs_quality"           ][:][:]; lwir_global_pcrs_quality=np.concatenate((lwir_global_pcrs_quality,source[imax,jmax]))
   source = lwva.variables["spatial_sample_quality"        ][:][:]; lwir_spatial_sample_quality=np.concatenate((lwir_spatial_sample_quality,source[imax,jmax]))
   source = mwva.variables["global_pc_scores"           ][:][:][:]; mwir_global_pc_scores=np.concatenate((mwir_global_pc_scores,source[imax,jmax]))
   source = lwva.variables["global_pc_scores"           ][:][:][:]; lwir_global_pc_scores=np.concatenate((lwir_global_pc_scores,source[imax,jmax]))

#---reconstruct radiances from this dwell

   rads=np.zeros((kmax,nchan))
   chan=np.zeros((kmax,nchan))

   for m in range(kmax):
      mm = kmax*(dwel-1)+m
      rads[m][lw0:lw1] = 100.* (np.dot(lwir_global_pc_scores[mm],recop_lw) + means_lw)
      rads[m][mw0:mw1] = 100.* (np.dot(mwir_global_pc_scores[mm],recop_mw) + means_mw)
      chan[m] = chans[:]

   radiance=np.concatenate((radiance,rads))
   chan_num=np.concatenate((chan_num,chan))

   ctime=0; ctime = tm.time(); ##print(btime-atime,ctime-btime); exit()

#---after all dwells are processed write into the container and dump the ioda file

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
container.add('chan_num', chan_num, ['*','*/RADCHN'])
container.add('radiance', radiance, ['*','*/RADCHN'])
netcdf.Encoder(description).encode(container,iodout)

