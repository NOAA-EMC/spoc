#!/usr/bin/env python3 # read mtg/irs netcdf files compute radiances and write to obsforge ioda file

import os
import sys
import bufr
import numpy as np
from netCDF4 import Dataset
from bufr.encoders import netcdf

# ----------------------------------------------------------------------------------------------------------
# define the apodising function
# ----------------------------------------------------------------------------------------------------------

def apply_hamming(x,nchan):
     hamming0 = np.float64(0.54); hamming1 = np.float64(0.23)
     y = x.copy(); ncn = nchan-1
     y[0] = (x[0]*hamming0+x[1]*hamming1)/(hamming0+hamming1)
     y[ncn] = (x[ncn]*hamming0+x[ncn-1]*hamming1)/(hamming0+hamming1)
     for icn in range(1,ncn):
          y[icn] = x[icn-1]*hamming1+x[icn]*hamming0+x[icn+1]*hamming1
     return y

# ----------------------------------------------------------------------------------------------------------
# read two arguments to define task
# ---------------------------------
# 1) the input data directory path
# 2) the output path/filename
# ----------------------------------------------------------------------------------------------------------

if len(sys.argv) < 3:
    print(f"{sys.argv[0]} needs <inpdir> and <iodout> ")
    exit()

inpdir=sys.argv[1]; print(sys.argv[1])
iodout=sys.argv[2]; print(sys.argv[2])

# ----------------------------------------------------------------------------------------------------------
# setup more filenames and parameters
# ----------------------------------------------------------------------------------------------------------

yamls   = "/scratch3/NCEPDEV/global/Jack.Woollen/spoc/dump/config/atmosphere"
parms   = "/scratch3/NCEPDEV/global/Jack.Woollen/spoc/dump/parm/atmosphere"
lwchans = os.path.join(parms,"irs_coopman_lw_channels.txt")
mwchans = os.path.join(parms,"irs_coopman_mw_channels.txt")
reconst = os.path.join(parms,"RSP_OPE_BASEEV_MTS1+IRS_20230925000000_V1_out.h5")
yaml    = os.path.join(yamls,"radiance_irs.yaml")

# ----------------------------------------------------------------------------------------------------------
# read the channel data for lw and mw ir wave numbers
# --------------------------------------------------

chan_lw = np.loadtxt(lwchans,dtype=int); chan_lw = chan_lw[chan_lw != 0]; chns_lw = chan_lw.shape[0]
chan_mw = np.loadtxt(mwchans,dtype=int); chan_mw = chan_mw[chan_mw != 0]; chns_mw = chan_mw.shape[0]

# ----------------------------------------------------------------------------------------------------------
# read reconstruction means and operators for lw and mw ir channels
# ----------------------------------------------------------------------------------------------------------

hd5file=Dataset(reconst)
means = hd5file['lwir'].variables['Mean'][:][:]; means_lw = means[chan_lw]
means = hd5file['mwir'].variables['Mean'][:][:]; means_mw = means[chan_mw]
recop = hd5file['lwir'].variables['ReconstructionOperator']; recop_lw = recop[:,chan_lw]; wnum_lw = recop.shape[1]
recop = hd5file['mwir'].variables['ReconstructionOperator']; recop_mw = recop[:,chan_mw]; wnum_mw = recop.shape[1]

# ----------------------------------------------------------------------------------------------------------
# hamming the eigenvectors
# ----------------------------------------------------------------------------------------------------------

means_lw = apply_hamming(means_lw,chns_lw)
means_mw = apply_hamming(means_mw,chns_mw)
for ipc in range(len(recop_lw)):
     recop_lw[ipc] = apply_hamming(recop_lw[ipc][:],chns_lw)
     recop_mw[ipc] = apply_hamming(recop_mw[ipc][:],chns_mw)

# ----------------------------------------------------------------------------------------------------------
# setup radiance channel output parameters 
# ----------------------------------------------------------------------------------------------------------

lw0 = 0
lw1 = chns_lw
mw0 = chns_lw
mw1 = mw0+chns_mw
nchan = chns_lw+chns_mw
chans = np.full(nchan,0,dtype=int)
chans[lw0:lw1] = chan_lw[:]
chans[mw0:mw1] = chan_mw[:]+wnum_lw

# ----------------------------------------------------------------------------------------------------------
# working parameters and arrays
# ----------------------------------------------------------------------------------------------------------

iter   = 1
itex   = -1         
imax   = np.zeros(1024,dtype=int)
jmax   = np.zeros(1024,dtype=int)
kmax   = 1024 
fill   = 1.e300

# ----------------------------------------------------------------------------------------------------------
# check that the files in the list exist and are readable
# ----------------------------------------------------------------------------------------------------------

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
   #print(iter, filename)
   try:
        irs  = Dataset(os.path.join(inpdir, filename))
        file_keep.append(filename)
   except OSError as e:
        sys.stderr.write(f"Skipping file '{filename}': failed to open as NetCDF: {e}\n")
        continue

# ----------------------------------------------------------------------------------------------------------
# create numpy array for all the dwell groups
# ----------------------------------------------------------------------------------------------------------

nfil=len(file_keep)      ; print('nfiles=',nfil)
nloc=len(file_keep)*kmax ; print('nlocs=',nloc)

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
chan_num=np.zeros((nloc,nchan))
radiance=np.zeros((nloc,nchan))

# ----------------------------------------------------------------------------------------------------------
# loop through the list of mtg-irs dwell files kept
# ----------------------------------------------------------------------------------------------------------

for filename in file_keep:
   print(iter,filename)
   irs  = Dataset(inpdir+"/"+filename)
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

# ----------------------------------------------------------------------------------------------------------
# select the hottest spot in each 5x5 box within the dwell
# ----------------------------------------------------------------------------------------------------------

   n = 0
   for a in range(0,160,5):
      for b in range(0,160,5):
         pcmax = -99.e99
         for i in range(1,4):
            for j in range(1,4):
               pc1=lwva.variables['global_pc_scores'][i+a][j+b][0]
               if abs(pc1) < fill: 
                  pcmax = max(pc1,pcmax)
                  if pcmax == pc1:
                     imax[n]=i+a
                     jmax[n]=j+b
         n=n+1

# ----------------------------------------------------------------------------------------------------------
# save the soundings selected from this dwell
# ----------------------------------------------------------------------------------------------------------

   m = (iter-1)*kmax-1
   for k in range(n):
      i = imax[k]
      j = jmax[k]
      m = m+1
      time[m]=loca.variables['time'][:]
      dwell_number[m] = loca.variables['dwell_number'][:]
      stroke_direction[m] = loca.variables['stroke_direction'][:]
      latitude[m] = loca.variables['latitude'][i][j]
      longitude[m] = loca.variables['longitude'][i][j]
      satellite_azimuth_angle[m] = loca.variables['satellite_azimuth_angle'][i][j]
      satellite_zenith_angle[m] = loca.variables['satellite_zenith_angle'][i][j]
      solar_azimuth_angle[m] = loca.variables['solar_azimuth_angle'][i][j]
      solar_zenith_angle[m] = loca.variables['solar_zenith_angle'][i][j]
      cloud_signal[m] = loca.variables['cloud_signal'][i][j]
      cloud_fraction[m] = loca.variables['cloud_fraction'][i][j]
      mwir_global_pc_scores[m] = mwva.variables['global_pc_scores'][i][j][:]
      mwir_global_pcr_scores[m] = mwva.variables['global_pcr_scores'][i][j]
      mwir_global_pcrs_quality[m] = mwva.variables['global_pcrs_quality'][i][j]
      mwir_spatial_sample_quality[m] = mwva.variables['spatial_sample_quality'][i][j]
      mwir_residual_energy[m] = mwva.variables['residual_energy'][:]
      lwir_global_pc_scores[m] = lwva.variables['global_pc_scores'][i][j][:]
      lwir_global_pcr_scores[m] = lwva.variables['global_pcr_scores'][i][j]
      lwir_global_pcrs_quality[m] = lwva.variables['global_pcrs_quality'][i][j]
      lwir_spatial_sample_quality[m] = lwva.variables['spatial_sample_quality'][i][j]
      lwir_residual_energy[m] = lwva.variables['residual_energy'][:]


      radiance[m][lw0:lw1] = 100.* (np.dot(lwir_global_pc_scores[m],recop_lw) + means_lw)
      radiance[m][mw0:mw1] = 100.* (np.dot(mwir_global_pc_scores[m],recop_mw) + means_mw)
      chan_num[m] = chans[:]

   if iter==itex:
      break
   iter=iter+1

#print(radiance.shape,lwir_global_pc_scores.shape,recop_lw.shape,means_lw.shape)
#print(radiance.dtype,lwir_global_pc_scores.dtype,recop_lw.dtype,means_lw.dtype)
#for i in range(300):
#   print(radiance[0][i])
#exit()

# ----------------------------------------------------------------------------------------------------------
# write into the container and the ioda dump file
# ----------------------------------------------------------------------------------------------------------

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

#for i in range(300):
#     print(i,radiance[0][i])

