#!/usr/bin/env python3 # read mtg/irs netcdf files compute radiances and write to obsforge ioda file

import os
import sys
import bufr
import time as tm
import numpy as np
import concurrent.futures
from netCDF4 import Dataset
from bufr.encoders import netcdf
from datetime import datetime, timezone
    
def main():

    # ---script to convert MTG-IRS data to IODA
    
    global file_keep,kmax,chunk,inpdir,thin,ind,jnd

    ATIME=0; BTIME=0; CTIME=0; DTIME=0; ETIME=0
    
    # ---read two arguments to define task
    # ---1) the input data directory path
    # ---2) the output path/filename
    # ---3) dir containing yamls        
    # ---4) dir continaing params    
    
    
    if len(sys.argv) < 4:
        print(f"{sys.argv[0]} needs <inpdir> <iodout> <yamls> <parms>")
        exit()
    
    inpdir = sys.argv[1]
    iodout = sys.argv[2]
    yamls  = sys.argv[3]
    parms  = sys.argv[4]
    
    print(inpdir)
    print(iodout)
    print(yamls)
    print(parms)
    
    #yaml = os.path.join(yamls, "pcscores_irs.yaml")
    yaml = os.path.join(yamls, "make-yaml/test.yaml")
    print(yaml)

    
    # ---check that the files in the list exist and are readable
    
    try:
        file_list = os.listdir(inpdir)
        file_keep = []
    except OSError as e:
        sys.stderr.write(f"Error reading input directory '{inpdir}': {e}n")
        sys.exit(1)
    
    netcdf_files = [f for f in file_list if f.lower().endswith(".nc")]
    if not netcdf_files:
        sys.stderr.write(f"No NetCDF files found in input directory '{inpdir}'.n")
        sys.exit(1)
    
    for filename in netcdf_files:
        try:
            #irs = Dataset(os.path.join(inpdir, filename))
            file_keep.append(filename)
            #break
            #irs.close()
        except OSError as e:
            sys.stderr.write(f"Skipping file '{filename}': failed to open as NetCDF: {e}n")
            continue
    
    dwel = 0
    nfil = len(file_keep)
    print('nfiles=', nfil)

    # ---select thinning to apply to  dwell files
    
    thin=5;  beg=1 ; end=4
    thin=10; beg=2 ; end=6
    thin=16; beg=4 ; end=10 
    thin=32; beg=12; end=18 
    
    gridSize = 160*160
    boxSize = thin*thin
    kmax = int(gridSize/boxSize)
    
    # ---create numpy arrays for selected channel data
    
    chan = 2000
    sensorCentralWavenumber = np.zeros(chan,dtype=np.float32)
    sensorChannelNumber     = np.zeros(chan,dtype=np.int32)
    
    # ---make a set of index offsets for the inner selection box
    
    n=0
    binsize=(end-beg)**2
    ind = np.zeros(binsize,int)
    jnd = np.zeros(binsize,int)
    for i in range(beg,end):
        for j in range(beg,end):
            ind[n] = i
            jnd[n] = j
            n=n+1
    
    # ---set up the file ranges for parallel processing
    
    nwork = 20
    nfils = len(file_keep)
    chunk = int(nfils/nwork) if int(nfils/nwork)*nwork == nfils else int(nfils/nwork)+1
    strt  = np.zeros(nwork,dtype=int)
    fini  = np.zeros(nwork,dtype=int)
    
    for idx in range(nwork):
        strt[idx] = idx*chunk
        fini[idx] = idx*chunk+chunk
    fini[nwork-1] = nfils
    
    # ---parallel loop through the list of mtg-irs dwell files to be processed
    
    rtime=tm.time()
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=nwork) as executor:
        data = np.array(list(executor.map(readloop,strt,fini)))
    
    # ---create numpy arrays to containerize the data in 
    
    vlist,dlist,nlist = varlist()
    
    # ---concatenate complete arrays from the parallel segments 
    
    mpic = 0
    for npic in range(nwork):
        kpic = np.count_nonzero(data[npic,0,:])
        vpic = 0
        for aray in range(len(vlist)):
            vlist[aray] = np.concatenate((vlist[aray],data[npic,vpic,0:kpic]))
            vpic = vpic+1
        mpic = mpic+kpic
    
    # ---make sure resultant arrays are the desired data types
    
    vpic = 0
    for aray in range(len(vlist)):
        vlist[aray] = vlist[aray].astype(dlist[vpic])
        vpic = vpic+1
    
    ATIME=tm.time()-rtime; rtime=tm.time()
    
    # ---containerize and encode the data  
    
    npic = 0
    container = bufr.DataContainer()
    description = bufr.encoders.Description(yaml)
    for aray in range(len(vlist)):
        container.add(nlist[npic], vlist[aray], ['*'])
        npic = npic + 1
    netcdf.Encoder(description).encode(container, iodout)
    
    BTIME=tm.time()-rtime
    
    # ---exit on time and in space
    print()
    print('thin/selec=',ATIME)
    print('ioda ncode=',BTIME)
    print('processing=',ATIME+BTIME+CTIME+DTIME+ETIME)
    print()
    print(mpic,' reports processed')

    return
    
# ---function to extract data from the mtg-irs dwell files to be processed

def readloop(beg,end):

    nloc = kmax*chunk 
    nvar = 26
    npcs = 150
    nwav = 2
    data = np.zeros((nvar+npcs*nwav,nloc),dtype=np.float32)
    datx = np.zeros((nwav,nloc,npcs),dtype=np.float32)
    mpic = 0

    #---process a range of input files
   
    print(beg,file_keep[beg])

    
    for fidx in range(beg,end):

        imax = np.full(kmax, -1, dtype=int)
        jmax = np.full(kmax, -1, dtype=int)
        kpic = 0
        n = 0

        try:
            irs = Dataset(os.path.join(inpdir, file_keep[fidx]))
        except OSError as e:
            break    

        loca = irs['data']
        mwva = irs['data/mwir/compressed']
        lwva = irs['data/lwir/compressed']

        #---process the qc criteria and zero pcscores in rejected spots

        pcsc = lwva.variables["global_pc_scores"][:][:][:]
        lons = loca.variables["longitude"][:][:]
        lats = loca.variables["latitude"][:][:]
        satz = loca.variables["satellite_zenith_angle"][:][:]

        fil = 2147483647
        qc1 = abs(pcsc[:,:,:]) <  fil
        lon = abs(lons) <= 180
        lat = abs(lats) <= 90 
        zan = abs(satz) <= 60 
        qc2 = lon & lat & zan

        pcsc = np.where (qc1, pcsc[:,:,:], 0)
        pcsc[:,:,0]  = np.where (qc2, pcsc[:,:,0], 0)

        #---select the hottest spot in each inner box in the dwell

        for a in range(0, 160, thin):
            for b in range(0, 160, thin):
                imx = ind+a
                jmx = jnd+b
                pc1 = np.abs(pcsc[imx,jmx,0])
                pcm = np.max(pc1)
                if pcm > 0:
                    inx = np.where(pc1 == pcm)[0]
                    imax[n] = imx[inx[0]]
                    jmax[n] = jmx[inx[0]]
                    n = n + 1
        if n == 0:
            continue 

        #---save the soundings selected from this dwell

        imax = imax[imax >= 0]
        jmax = jmax[jmax >= 0]
        kpic = imax.shape[0]
        lpic = mpic        
        mpic = mpic+kpic       

        scalar = np.zeros(kpic)

        scalar = loca.variables["time"][:]   
        data[0,lpic:mpic] = scalar

        scalar = loca.variables["dwell_number"][:]   
        data[1,lpic:mpic] = scalar

        scalar = loca.variables["dwell_type"][:]   
        data[2,lpic:mpic] = scalar

        source = lwva.variables["detector_sample_quality"][:][:]
        data[3,lpic:mpic] = source[imax,jmax]

        source = mwva.variables["detector_sample_quality"][:][:]
        data[4,lpic:mpic] = source[imax,jmax]

        source = lwva.variables["global_pcrs_quality"][:][:]
        data[5,lpic:mpic] = source[imax,jmax]

        source = mwva.variables["global_pcrs_quality"][:][:]
        data[6,lpic:mpic] = source[imax,jmax]

        source = lwva.variables["spatial_sample_quality"][:][:]
        data[7,lpic:mpic] = source[imax,jmax]

        source = mwva.variables["spatial_sample_quality"][:][:]
        data[8,lpic:mpic] = source[imax,jmax]

        data[9,lpic:mpic] = 0     

        data[10,lpic:mpic] = 0     

        data[11,lpic:mpic] = 0     

        data[12,lpic:mpic] = 0     

        data[13,lpic:mpic] = 0     

        source = loca.variables["cloud_fraction"][:][:]
        data[14,lpic:mpic] = source[imax,jmax]

        source = loca.variables["cloud_signal"][:][:]
        data[15,lpic:mpic] = source[imax,jmax]

        source = loca.variables["dust_warning"][:][:]
        data[16,lpic:mpic] = source[imax,jmax]

        source = loca.variables["latitude"][:][:]
        data[17,lpic:mpic] = source[imax,jmax]

        source = loca.variables["longitude"][:][:]
        data[18,lpic:mpic] = source[imax,jmax]

        source = loca.variables["satellite_azimuth_angle"][:][:]
        data[19,lpic:mpic] = source[imax,jmax]

        source = loca.variables["satellite_zenith_angle"][:][:]
        data[20,lpic:mpic] = source[imax,jmax]

        source = loca.variables["solar_azimuth_angle"][:][:]
        data[21,lpic:mpic] = source[imax,jmax]

        source = loca.variables["solar_zenith_angle"][:][:]
        data[22,lpic:mpic] = source[imax,jmax]

        data[23,lpic:mpic] = 0     

        source = lwva.variables["global_pcr_scores"][:][:]
        data[24,lpic:mpic] = source[imax,jmax]

        source = mwva.variables["global_pcr_scores"][:][:]
        data[25,lpic:mpic] = source[imax,jmax]

        source = lwva.variables["global_pc_scores"][:][:][:]
        datx[0,lpic:mpic] = source[imax,jmax]

        source = mwva.variables["global_pc_scores"][:][:][:]
        datx[1,lpic:mpic] = source[imax,jmax]

        for ipc in range(npcs):
            lwpc = nvar + ipc
            mwpc = lwpc + npcs
            data[lwpc,:] = datx[0,:,ipc]
            data[mwpc,:] = datx[1,:,ipc]

        irs.close()

    return data
  
def varlist():
    # 1. Define specific metadata variables
    # Format: (name, dtype)
    meta_config = [
        ('dateTime', np.int64),
        ('dwellNumber', np.int32),
        ('dwellType', np.int32),
        ('detectorSampleQualityLw', np.int32),
        ('detectorSampleQualityMw', np.int32),
        ('globalPcrQualityLw', np.int32),
        ('globalPcrQualityMw', np.int32),
        ('spatialSampleQualityLw', np.int32),
        ('spatialSampleQualityMw', np.int32),
        ('overallQuality', np.int32),
        ('satelliteIdentifier', np.int32),
        ('sensorScanPosition', np.int32),
        ('sensorChannelNumber', np.int32),
        ('sensorCentralWavenumber', np.float32),
        ('cloudFraction', np.float32),
        ('cloudSignal', np.float32),
        ('dustWarning', np.float32),
        ('latitude', np.float32),
        ('longitude', np.float32),
        ('sensorAzimuthAngle', np.float32),
        ('sensorZenithAngle', np.float32),
        ('solarAzimuthAngle', np.float32),
        ('solarZenithAngle', np.float32),
        ('sensorViewAngle', np.float32),
        ('globalPcrScoresLw', np.float32),
        ('globalPcrScoresMw', np.float32),
    ]

    vlist, dlist, nlist = [], [], []

    # 2. Add the metadata variables to the lists
    for name, dtype in meta_config:
        nlist.append(name)
        dlist.append(dtype)
        vlist.append(np.zeros(0, dtype=dtype))

    # 3. Programmatically add the 300 Principal Component Scores
    for i in range(1, 301):
        name = f'principalComponentScore{i}'
        nlist.append(name)
        dlist.append(np.float32)
        vlist.append(np.zeros(0, dtype=np.float32))

    return vlist, dlist, nlist

if __name__ == "__main__":
    main()

