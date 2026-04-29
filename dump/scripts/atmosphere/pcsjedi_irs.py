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
    global vlist,dlist,nlist
    vlist = []
    dlist = []
    nlist = []
    global dateTime;dateTime=np.zeros(0);vlist.append(dateTime);dlist.append(np.int64);nlist.append('dateTime')
    global dwellNumber;dwellNumber=np.zeros(0);vlist.append(dwellNumber);dlist.append(np.int32);nlist.append('dwellNumber')
    global dwellType;dwellType=np.zeros(0);vlist.append(dwellType);dlist.append(np.int32);nlist.append('dwellType')
    global detectorSampleQualityLw;detectorSampleQualityLw=np.zeros(0);vlist.append(detectorSampleQualityLw);dlist.append(np.int32);nlist.append('detectorSampleQualityLw')
    global detectorSampleQualityMw;detectorSampleQualityMw=np.zeros(0);vlist.append(detectorSampleQualityMw);dlist.append(np.int32);nlist.append('detectorSampleQualityMw')
    global globalPcrQualityLw;globalPcrQualityLw=np.zeros(0);vlist.append(globalPcrQualityLw);dlist.append(np.int32);nlist.append('globalPcrQualityLw')
    global globalPcrQualityMw;globalPcrQualityMw=np.zeros(0);vlist.append(globalPcrQualityMw);dlist.append(np.int32);nlist.append('globalPcrQualityMw')
    global spatialSampleQualityLw;spatialSampleQualityLw=np.zeros(0);vlist.append(spatialSampleQualityLw);dlist.append(np.int32);nlist.append('spatialSampleQualityLw')
    global spatialSampleQualityMw;spatialSampleQualityMw=np.zeros(0);vlist.append(spatialSampleQualityMw);dlist.append(np.int32);nlist.append('spatialSampleQualityMw')
    global overallQuality;overallQuality=np.zeros(0);vlist.append(overallQuality);dlist.append(np.int32);nlist.append('overallQuality')
    global satelliteIdentifier;satelliteIdentifier=np.zeros(0);vlist.append(satelliteIdentifier);dlist.append(np.int32);nlist.append('satelliteIdentifier')
    global sensorScanPosition;sensorScanPosition=np.zeros(0);vlist.append(sensorScanPosition);dlist.append(np.int32);nlist.append('sensorScanPosition')
    global sensorChannelNumber;sensorChannelNumber=np.zeros(0);vlist.append(sensorChannelNumber);dlist.append(np.int32);nlist.append('sensorChannelNumber')
    global sensorCentralWavenumber;sensorCentralWavenumber=np.zeros(0);vlist.append(sensorCentralWavenumber);dlist.append(np.float32);nlist.append('sensorCentralWavenumber')
    global cloudFraction;cloudFraction=np.zeros(0);vlist.append(cloudFraction);dlist.append(np.float32);nlist.append('cloudFraction')
    global cloudSignal;cloudSignal=np.zeros(0);vlist.append(cloudSignal);dlist.append(np.float32);nlist.append('cloudSignal')
    global dustWarning;dustWarning=np.zeros(0);vlist.append(dustWarning);dlist.append(np.float32);nlist.append('dustWarning')
    global latitude;latitude=np.zeros(0);vlist.append(latitude);dlist.append(np.float32);nlist.append('latitude')
    global longitude;longitude=np.zeros(0);vlist.append(longitude);dlist.append(np.float32);nlist.append('longitude')
    global sensorAzimuthAngle;sensorAzimuthAngle=np.zeros(0);vlist.append(sensorAzimuthAngle);dlist.append(np.float32);nlist.append('sensorAzimuthAngle')
    global sensorZenithAngle;sensorZenithAngle=np.zeros(0);vlist.append(sensorZenithAngle);dlist.append(np.float32);nlist.append('sensorZenithAngle')
    global solarAzimuthAngle;solarAzimuthAngle=np.zeros(0);vlist.append(solarAzimuthAngle);dlist.append(np.float32);nlist.append('solarAzimuthAngle')
    global solarZenithAngle;solarZenithAngle=np.zeros(0);vlist.append(solarZenithAngle);dlist.append(np.float32);nlist.append('solarZenithAngle')
    global sensorViewAngle;sensorViewAngle=np.zeros(0);vlist.append(sensorViewAngle);dlist.append(np.float32);nlist.append('sensorViewAngle')
    global globalPcrScoresLw;globalPcrScoresLw=np.zeros(0);vlist.append(globalPcrScoresLw);dlist.append(np.float32);nlist.append('globalPcrScoresLw')
    global globalPcrScoresMw;globalPcrScoresMw=np.zeros(0);vlist.append(globalPcrScoresMw);dlist.append(np.float32);nlist.append('globalPcrScoresMw')
    global principalComponentScore1;principalComponentScore1=np.zeros(0);vlist.append(principalComponentScore1);dlist.append(np.float32);nlist.append('principalComponentScore1')
    global principalComponentScore2;principalComponentScore2=np.zeros(0);vlist.append(principalComponentScore2);dlist.append(np.float32);nlist.append('principalComponentScore2')
    global principalComponentScore3;principalComponentScore3=np.zeros(0);vlist.append(principalComponentScore3);dlist.append(np.float32);nlist.append('principalComponentScore3')
    global principalComponentScore4;principalComponentScore4=np.zeros(0);vlist.append(principalComponentScore4);dlist.append(np.float32);nlist.append('principalComponentScore4')
    global principalComponentScore5;principalComponentScore5=np.zeros(0);vlist.append(principalComponentScore5);dlist.append(np.float32);nlist.append('principalComponentScore5')
    global principalComponentScore6;principalComponentScore6=np.zeros(0);vlist.append(principalComponentScore6);dlist.append(np.float32);nlist.append('principalComponentScore6')
    global principalComponentScore7;principalComponentScore7=np.zeros(0);vlist.append(principalComponentScore7);dlist.append(np.float32);nlist.append('principalComponentScore7')
    global principalComponentScore8;principalComponentScore8=np.zeros(0);vlist.append(principalComponentScore8);dlist.append(np.float32);nlist.append('principalComponentScore8')
    global principalComponentScore9;principalComponentScore9=np.zeros(0);vlist.append(principalComponentScore9);dlist.append(np.float32);nlist.append('principalComponentScore9')
    global principalComponentScore10;principalComponentScore10=np.zeros(0);vlist.append(principalComponentScore10);dlist.append(np.float32);nlist.append('principalComponentScore10')
    global principalComponentScore11;principalComponentScore11=np.zeros(0);vlist.append(principalComponentScore11);dlist.append(np.float32);nlist.append('principalComponentScore11')
    global principalComponentScore12;principalComponentScore12=np.zeros(0);vlist.append(principalComponentScore12);dlist.append(np.float32);nlist.append('principalComponentScore12')
    global principalComponentScore13;principalComponentScore13=np.zeros(0);vlist.append(principalComponentScore13);dlist.append(np.float32);nlist.append('principalComponentScore13')
    global principalComponentScore14;principalComponentScore14=np.zeros(0);vlist.append(principalComponentScore14);dlist.append(np.float32);nlist.append('principalComponentScore14')
    global principalComponentScore15;principalComponentScore15=np.zeros(0);vlist.append(principalComponentScore15);dlist.append(np.float32);nlist.append('principalComponentScore15')
    global principalComponentScore16;principalComponentScore16=np.zeros(0);vlist.append(principalComponentScore16);dlist.append(np.float32);nlist.append('principalComponentScore16')
    global principalComponentScore17;principalComponentScore17=np.zeros(0);vlist.append(principalComponentScore17);dlist.append(np.float32);nlist.append('principalComponentScore17')
    global principalComponentScore18;principalComponentScore18=np.zeros(0);vlist.append(principalComponentScore18);dlist.append(np.float32);nlist.append('principalComponentScore18')
    global principalComponentScore19;principalComponentScore19=np.zeros(0);vlist.append(principalComponentScore19);dlist.append(np.float32);nlist.append('principalComponentScore19')
    global principalComponentScore20;principalComponentScore20=np.zeros(0);vlist.append(principalComponentScore20);dlist.append(np.float32);nlist.append('principalComponentScore20')
    global principalComponentScore21;principalComponentScore21=np.zeros(0);vlist.append(principalComponentScore21);dlist.append(np.float32);nlist.append('principalComponentScore21')
    global principalComponentScore22;principalComponentScore22=np.zeros(0);vlist.append(principalComponentScore22);dlist.append(np.float32);nlist.append('principalComponentScore22')
    global principalComponentScore23;principalComponentScore23=np.zeros(0);vlist.append(principalComponentScore23);dlist.append(np.float32);nlist.append('principalComponentScore23')
    global principalComponentScore24;principalComponentScore24=np.zeros(0);vlist.append(principalComponentScore24);dlist.append(np.float32);nlist.append('principalComponentScore24')
    global principalComponentScore25;principalComponentScore25=np.zeros(0);vlist.append(principalComponentScore25);dlist.append(np.float32);nlist.append('principalComponentScore25')
    global principalComponentScore26;principalComponentScore26=np.zeros(0);vlist.append(principalComponentScore26);dlist.append(np.float32);nlist.append('principalComponentScore26')
    global principalComponentScore27;principalComponentScore27=np.zeros(0);vlist.append(principalComponentScore27);dlist.append(np.float32);nlist.append('principalComponentScore27')
    global principalComponentScore28;principalComponentScore28=np.zeros(0);vlist.append(principalComponentScore28);dlist.append(np.float32);nlist.append('principalComponentScore28')
    global principalComponentScore29;principalComponentScore29=np.zeros(0);vlist.append(principalComponentScore29);dlist.append(np.float32);nlist.append('principalComponentScore29')
    global principalComponentScore30;principalComponentScore30=np.zeros(0);vlist.append(principalComponentScore30);dlist.append(np.float32);nlist.append('principalComponentScore30')
    global principalComponentScore31;principalComponentScore31=np.zeros(0);vlist.append(principalComponentScore31);dlist.append(np.float32);nlist.append('principalComponentScore31')
    global principalComponentScore32;principalComponentScore32=np.zeros(0);vlist.append(principalComponentScore32);dlist.append(np.float32);nlist.append('principalComponentScore32')
    global principalComponentScore33;principalComponentScore33=np.zeros(0);vlist.append(principalComponentScore33);dlist.append(np.float32);nlist.append('principalComponentScore33')
    global principalComponentScore34;principalComponentScore34=np.zeros(0);vlist.append(principalComponentScore34);dlist.append(np.float32);nlist.append('principalComponentScore34')
    global principalComponentScore35;principalComponentScore35=np.zeros(0);vlist.append(principalComponentScore35);dlist.append(np.float32);nlist.append('principalComponentScore35')
    global principalComponentScore36;principalComponentScore36=np.zeros(0);vlist.append(principalComponentScore36);dlist.append(np.float32);nlist.append('principalComponentScore36')
    global principalComponentScore37;principalComponentScore37=np.zeros(0);vlist.append(principalComponentScore37);dlist.append(np.float32);nlist.append('principalComponentScore37')
    global principalComponentScore38;principalComponentScore38=np.zeros(0);vlist.append(principalComponentScore38);dlist.append(np.float32);nlist.append('principalComponentScore38')
    global principalComponentScore39;principalComponentScore39=np.zeros(0);vlist.append(principalComponentScore39);dlist.append(np.float32);nlist.append('principalComponentScore39')
    global principalComponentScore40;principalComponentScore40=np.zeros(0);vlist.append(principalComponentScore40);dlist.append(np.float32);nlist.append('principalComponentScore40')
    global principalComponentScore41;principalComponentScore41=np.zeros(0);vlist.append(principalComponentScore41);dlist.append(np.float32);nlist.append('principalComponentScore41')
    global principalComponentScore42;principalComponentScore42=np.zeros(0);vlist.append(principalComponentScore42);dlist.append(np.float32);nlist.append('principalComponentScore42')
    global principalComponentScore43;principalComponentScore43=np.zeros(0);vlist.append(principalComponentScore43);dlist.append(np.float32);nlist.append('principalComponentScore43')
    global principalComponentScore44;principalComponentScore44=np.zeros(0);vlist.append(principalComponentScore44);dlist.append(np.float32);nlist.append('principalComponentScore44')
    global principalComponentScore45;principalComponentScore45=np.zeros(0);vlist.append(principalComponentScore45);dlist.append(np.float32);nlist.append('principalComponentScore45')
    global principalComponentScore46;principalComponentScore46=np.zeros(0);vlist.append(principalComponentScore46);dlist.append(np.float32);nlist.append('principalComponentScore46')
    global principalComponentScore47;principalComponentScore47=np.zeros(0);vlist.append(principalComponentScore47);dlist.append(np.float32);nlist.append('principalComponentScore47')
    global principalComponentScore48;principalComponentScore48=np.zeros(0);vlist.append(principalComponentScore48);dlist.append(np.float32);nlist.append('principalComponentScore48')
    global principalComponentScore49;principalComponentScore49=np.zeros(0);vlist.append(principalComponentScore49);dlist.append(np.float32);nlist.append('principalComponentScore49')
    global principalComponentScore50;principalComponentScore50=np.zeros(0);vlist.append(principalComponentScore50);dlist.append(np.float32);nlist.append('principalComponentScore50')
    global principalComponentScore51;principalComponentScore51=np.zeros(0);vlist.append(principalComponentScore51);dlist.append(np.float32);nlist.append('principalComponentScore51')
    global principalComponentScore52;principalComponentScore52=np.zeros(0);vlist.append(principalComponentScore52);dlist.append(np.float32);nlist.append('principalComponentScore52')
    global principalComponentScore53;principalComponentScore53=np.zeros(0);vlist.append(principalComponentScore53);dlist.append(np.float32);nlist.append('principalComponentScore53')
    global principalComponentScore54;principalComponentScore54=np.zeros(0);vlist.append(principalComponentScore54);dlist.append(np.float32);nlist.append('principalComponentScore54')
    global principalComponentScore55;principalComponentScore55=np.zeros(0);vlist.append(principalComponentScore55);dlist.append(np.float32);nlist.append('principalComponentScore55')
    global principalComponentScore56;principalComponentScore56=np.zeros(0);vlist.append(principalComponentScore56);dlist.append(np.float32);nlist.append('principalComponentScore56')
    global principalComponentScore57;principalComponentScore57=np.zeros(0);vlist.append(principalComponentScore57);dlist.append(np.float32);nlist.append('principalComponentScore57')
    global principalComponentScore58;principalComponentScore58=np.zeros(0);vlist.append(principalComponentScore58);dlist.append(np.float32);nlist.append('principalComponentScore58')
    global principalComponentScore59;principalComponentScore59=np.zeros(0);vlist.append(principalComponentScore59);dlist.append(np.float32);nlist.append('principalComponentScore59')
    global principalComponentScore60;principalComponentScore60=np.zeros(0);vlist.append(principalComponentScore60);dlist.append(np.float32);nlist.append('principalComponentScore60')
    global principalComponentScore61;principalComponentScore61=np.zeros(0);vlist.append(principalComponentScore61);dlist.append(np.float32);nlist.append('principalComponentScore61')
    global principalComponentScore62;principalComponentScore62=np.zeros(0);vlist.append(principalComponentScore62);dlist.append(np.float32);nlist.append('principalComponentScore62')
    global principalComponentScore63;principalComponentScore63=np.zeros(0);vlist.append(principalComponentScore63);dlist.append(np.float32);nlist.append('principalComponentScore63')
    global principalComponentScore64;principalComponentScore64=np.zeros(0);vlist.append(principalComponentScore64);dlist.append(np.float32);nlist.append('principalComponentScore64')
    global principalComponentScore65;principalComponentScore65=np.zeros(0);vlist.append(principalComponentScore65);dlist.append(np.float32);nlist.append('principalComponentScore65')
    global principalComponentScore66;principalComponentScore66=np.zeros(0);vlist.append(principalComponentScore66);dlist.append(np.float32);nlist.append('principalComponentScore66')
    global principalComponentScore67;principalComponentScore67=np.zeros(0);vlist.append(principalComponentScore67);dlist.append(np.float32);nlist.append('principalComponentScore67')
    global principalComponentScore68;principalComponentScore68=np.zeros(0);vlist.append(principalComponentScore68);dlist.append(np.float32);nlist.append('principalComponentScore68')
    global principalComponentScore69;principalComponentScore69=np.zeros(0);vlist.append(principalComponentScore69);dlist.append(np.float32);nlist.append('principalComponentScore69')
    global principalComponentScore70;principalComponentScore70=np.zeros(0);vlist.append(principalComponentScore70);dlist.append(np.float32);nlist.append('principalComponentScore70')
    global principalComponentScore71;principalComponentScore71=np.zeros(0);vlist.append(principalComponentScore71);dlist.append(np.float32);nlist.append('principalComponentScore71')
    global principalComponentScore72;principalComponentScore72=np.zeros(0);vlist.append(principalComponentScore72);dlist.append(np.float32);nlist.append('principalComponentScore72')
    global principalComponentScore73;principalComponentScore73=np.zeros(0);vlist.append(principalComponentScore73);dlist.append(np.float32);nlist.append('principalComponentScore73')
    global principalComponentScore74;principalComponentScore74=np.zeros(0);vlist.append(principalComponentScore74);dlist.append(np.float32);nlist.append('principalComponentScore74')
    global principalComponentScore75;principalComponentScore75=np.zeros(0);vlist.append(principalComponentScore75);dlist.append(np.float32);nlist.append('principalComponentScore75')
    global principalComponentScore76;principalComponentScore76=np.zeros(0);vlist.append(principalComponentScore76);dlist.append(np.float32);nlist.append('principalComponentScore76')
    global principalComponentScore77;principalComponentScore77=np.zeros(0);vlist.append(principalComponentScore77);dlist.append(np.float32);nlist.append('principalComponentScore77')
    global principalComponentScore78;principalComponentScore78=np.zeros(0);vlist.append(principalComponentScore78);dlist.append(np.float32);nlist.append('principalComponentScore78')
    global principalComponentScore79;principalComponentScore79=np.zeros(0);vlist.append(principalComponentScore79);dlist.append(np.float32);nlist.append('principalComponentScore79')
    global principalComponentScore80;principalComponentScore80=np.zeros(0);vlist.append(principalComponentScore80);dlist.append(np.float32);nlist.append('principalComponentScore80')
    global principalComponentScore81;principalComponentScore81=np.zeros(0);vlist.append(principalComponentScore81);dlist.append(np.float32);nlist.append('principalComponentScore81')
    global principalComponentScore82;principalComponentScore82=np.zeros(0);vlist.append(principalComponentScore82);dlist.append(np.float32);nlist.append('principalComponentScore82')
    global principalComponentScore83;principalComponentScore83=np.zeros(0);vlist.append(principalComponentScore83);dlist.append(np.float32);nlist.append('principalComponentScore83')
    global principalComponentScore84;principalComponentScore84=np.zeros(0);vlist.append(principalComponentScore84);dlist.append(np.float32);nlist.append('principalComponentScore84')
    global principalComponentScore85;principalComponentScore85=np.zeros(0);vlist.append(principalComponentScore85);dlist.append(np.float32);nlist.append('principalComponentScore85')
    global principalComponentScore86;principalComponentScore86=np.zeros(0);vlist.append(principalComponentScore86);dlist.append(np.float32);nlist.append('principalComponentScore86')
    global principalComponentScore87;principalComponentScore87=np.zeros(0);vlist.append(principalComponentScore87);dlist.append(np.float32);nlist.append('principalComponentScore87')
    global principalComponentScore88;principalComponentScore88=np.zeros(0);vlist.append(principalComponentScore88);dlist.append(np.float32);nlist.append('principalComponentScore88')
    global principalComponentScore89;principalComponentScore89=np.zeros(0);vlist.append(principalComponentScore89);dlist.append(np.float32);nlist.append('principalComponentScore89')
    global principalComponentScore90;principalComponentScore90=np.zeros(0);vlist.append(principalComponentScore90);dlist.append(np.float32);nlist.append('principalComponentScore90')
    global principalComponentScore91;principalComponentScore91=np.zeros(0);vlist.append(principalComponentScore91);dlist.append(np.float32);nlist.append('principalComponentScore91')
    global principalComponentScore92;principalComponentScore92=np.zeros(0);vlist.append(principalComponentScore92);dlist.append(np.float32);nlist.append('principalComponentScore92')
    global principalComponentScore93;principalComponentScore93=np.zeros(0);vlist.append(principalComponentScore93);dlist.append(np.float32);nlist.append('principalComponentScore93')
    global principalComponentScore94;principalComponentScore94=np.zeros(0);vlist.append(principalComponentScore94);dlist.append(np.float32);nlist.append('principalComponentScore94')
    global principalComponentScore95;principalComponentScore95=np.zeros(0);vlist.append(principalComponentScore95);dlist.append(np.float32);nlist.append('principalComponentScore95')
    global principalComponentScore96;principalComponentScore96=np.zeros(0);vlist.append(principalComponentScore96);dlist.append(np.float32);nlist.append('principalComponentScore96')
    global principalComponentScore97;principalComponentScore97=np.zeros(0);vlist.append(principalComponentScore97);dlist.append(np.float32);nlist.append('principalComponentScore97')
    global principalComponentScore98;principalComponentScore98=np.zeros(0);vlist.append(principalComponentScore98);dlist.append(np.float32);nlist.append('principalComponentScore98')
    global principalComponentScore99;principalComponentScore99=np.zeros(0);vlist.append(principalComponentScore99);dlist.append(np.float32);nlist.append('principalComponentScore99')
    global principalComponentScore100;principalComponentScore100=np.zeros(0);vlist.append(principalComponentScore100);dlist.append(np.float32);nlist.append('principalComponentScore100')
    global principalComponentScore101;principalComponentScore101=np.zeros(0);vlist.append(principalComponentScore101);dlist.append(np.float32);nlist.append('principalComponentScore101')
    global principalComponentScore102;principalComponentScore102=np.zeros(0);vlist.append(principalComponentScore102);dlist.append(np.float32);nlist.append('principalComponentScore102')
    global principalComponentScore103;principalComponentScore103=np.zeros(0);vlist.append(principalComponentScore103);dlist.append(np.float32);nlist.append('principalComponentScore103')
    global principalComponentScore104;principalComponentScore104=np.zeros(0);vlist.append(principalComponentScore104);dlist.append(np.float32);nlist.append('principalComponentScore104')
    global principalComponentScore105;principalComponentScore105=np.zeros(0);vlist.append(principalComponentScore105);dlist.append(np.float32);nlist.append('principalComponentScore105')
    global principalComponentScore106;principalComponentScore106=np.zeros(0);vlist.append(principalComponentScore106);dlist.append(np.float32);nlist.append('principalComponentScore106')
    global principalComponentScore107;principalComponentScore107=np.zeros(0);vlist.append(principalComponentScore107);dlist.append(np.float32);nlist.append('principalComponentScore107')
    global principalComponentScore108;principalComponentScore108=np.zeros(0);vlist.append(principalComponentScore108);dlist.append(np.float32);nlist.append('principalComponentScore108')
    global principalComponentScore109;principalComponentScore109=np.zeros(0);vlist.append(principalComponentScore109);dlist.append(np.float32);nlist.append('principalComponentScore109')
    global principalComponentScore110;principalComponentScore110=np.zeros(0);vlist.append(principalComponentScore110);dlist.append(np.float32);nlist.append('principalComponentScore110')
    global principalComponentScore111;principalComponentScore111=np.zeros(0);vlist.append(principalComponentScore111);dlist.append(np.float32);nlist.append('principalComponentScore111')
    global principalComponentScore112;principalComponentScore112=np.zeros(0);vlist.append(principalComponentScore112);dlist.append(np.float32);nlist.append('principalComponentScore112')
    global principalComponentScore113;principalComponentScore113=np.zeros(0);vlist.append(principalComponentScore113);dlist.append(np.float32);nlist.append('principalComponentScore113')
    global principalComponentScore114;principalComponentScore114=np.zeros(0);vlist.append(principalComponentScore114);dlist.append(np.float32);nlist.append('principalComponentScore114')
    global principalComponentScore115;principalComponentScore115=np.zeros(0);vlist.append(principalComponentScore115);dlist.append(np.float32);nlist.append('principalComponentScore115')
    global principalComponentScore116;principalComponentScore116=np.zeros(0);vlist.append(principalComponentScore116);dlist.append(np.float32);nlist.append('principalComponentScore116')
    global principalComponentScore117;principalComponentScore117=np.zeros(0);vlist.append(principalComponentScore117);dlist.append(np.float32);nlist.append('principalComponentScore117')
    global principalComponentScore118;principalComponentScore118=np.zeros(0);vlist.append(principalComponentScore118);dlist.append(np.float32);nlist.append('principalComponentScore118')
    global principalComponentScore119;principalComponentScore119=np.zeros(0);vlist.append(principalComponentScore119);dlist.append(np.float32);nlist.append('principalComponentScore119')
    global principalComponentScore120;principalComponentScore120=np.zeros(0);vlist.append(principalComponentScore120);dlist.append(np.float32);nlist.append('principalComponentScore120')
    global principalComponentScore121;principalComponentScore121=np.zeros(0);vlist.append(principalComponentScore121);dlist.append(np.float32);nlist.append('principalComponentScore121')
    global principalComponentScore122;principalComponentScore122=np.zeros(0);vlist.append(principalComponentScore122);dlist.append(np.float32);nlist.append('principalComponentScore122')
    global principalComponentScore123;principalComponentScore123=np.zeros(0);vlist.append(principalComponentScore123);dlist.append(np.float32);nlist.append('principalComponentScore123')
    global principalComponentScore124;principalComponentScore124=np.zeros(0);vlist.append(principalComponentScore124);dlist.append(np.float32);nlist.append('principalComponentScore124')
    global principalComponentScore125;principalComponentScore125=np.zeros(0);vlist.append(principalComponentScore125);dlist.append(np.float32);nlist.append('principalComponentScore125')
    global principalComponentScore126;principalComponentScore126=np.zeros(0);vlist.append(principalComponentScore126);dlist.append(np.float32);nlist.append('principalComponentScore126')
    global principalComponentScore127;principalComponentScore127=np.zeros(0);vlist.append(principalComponentScore127);dlist.append(np.float32);nlist.append('principalComponentScore127')
    global principalComponentScore128;principalComponentScore128=np.zeros(0);vlist.append(principalComponentScore128);dlist.append(np.float32);nlist.append('principalComponentScore128')
    global principalComponentScore129;principalComponentScore129=np.zeros(0);vlist.append(principalComponentScore129);dlist.append(np.float32);nlist.append('principalComponentScore129')
    global principalComponentScore130;principalComponentScore130=np.zeros(0);vlist.append(principalComponentScore130);dlist.append(np.float32);nlist.append('principalComponentScore130')
    global principalComponentScore131;principalComponentScore131=np.zeros(0);vlist.append(principalComponentScore131);dlist.append(np.float32);nlist.append('principalComponentScore131')
    global principalComponentScore132;principalComponentScore132=np.zeros(0);vlist.append(principalComponentScore132);dlist.append(np.float32);nlist.append('principalComponentScore132')
    global principalComponentScore133;principalComponentScore133=np.zeros(0);vlist.append(principalComponentScore133);dlist.append(np.float32);nlist.append('principalComponentScore133')
    global principalComponentScore134;principalComponentScore134=np.zeros(0);vlist.append(principalComponentScore134);dlist.append(np.float32);nlist.append('principalComponentScore134')
    global principalComponentScore135;principalComponentScore135=np.zeros(0);vlist.append(principalComponentScore135);dlist.append(np.float32);nlist.append('principalComponentScore135')
    global principalComponentScore136;principalComponentScore136=np.zeros(0);vlist.append(principalComponentScore136);dlist.append(np.float32);nlist.append('principalComponentScore136')
    global principalComponentScore137;principalComponentScore137=np.zeros(0);vlist.append(principalComponentScore137);dlist.append(np.float32);nlist.append('principalComponentScore137')
    global principalComponentScore138;principalComponentScore138=np.zeros(0);vlist.append(principalComponentScore138);dlist.append(np.float32);nlist.append('principalComponentScore138')
    global principalComponentScore139;principalComponentScore139=np.zeros(0);vlist.append(principalComponentScore139);dlist.append(np.float32);nlist.append('principalComponentScore139')
    global principalComponentScore140;principalComponentScore140=np.zeros(0);vlist.append(principalComponentScore140);dlist.append(np.float32);nlist.append('principalComponentScore140')
    global principalComponentScore141;principalComponentScore141=np.zeros(0);vlist.append(principalComponentScore141);dlist.append(np.float32);nlist.append('principalComponentScore141')
    global principalComponentScore142;principalComponentScore142=np.zeros(0);vlist.append(principalComponentScore142);dlist.append(np.float32);nlist.append('principalComponentScore142')
    global principalComponentScore143;principalComponentScore143=np.zeros(0);vlist.append(principalComponentScore143);dlist.append(np.float32);nlist.append('principalComponentScore143')
    global principalComponentScore144;principalComponentScore144=np.zeros(0);vlist.append(principalComponentScore144);dlist.append(np.float32);nlist.append('principalComponentScore144')
    global principalComponentScore145;principalComponentScore145=np.zeros(0);vlist.append(principalComponentScore145);dlist.append(np.float32);nlist.append('principalComponentScore145')
    global principalComponentScore146;principalComponentScore146=np.zeros(0);vlist.append(principalComponentScore146);dlist.append(np.float32);nlist.append('principalComponentScore146')
    global principalComponentScore147;principalComponentScore147=np.zeros(0);vlist.append(principalComponentScore147);dlist.append(np.float32);nlist.append('principalComponentScore147')
    global principalComponentScore148;principalComponentScore148=np.zeros(0);vlist.append(principalComponentScore148);dlist.append(np.float32);nlist.append('principalComponentScore148')
    global principalComponentScore149;principalComponentScore149=np.zeros(0);vlist.append(principalComponentScore149);dlist.append(np.float32);nlist.append('principalComponentScore149')
    global principalComponentScore150;principalComponentScore150=np.zeros(0);vlist.append(principalComponentScore150);dlist.append(np.float32);nlist.append('principalComponentScore150')
    global principalComponentScore151;principalComponentScore151=np.zeros(0);vlist.append(principalComponentScore151);dlist.append(np.float32);nlist.append('principalComponentScore151')
    global principalComponentScore152;principalComponentScore152=np.zeros(0);vlist.append(principalComponentScore152);dlist.append(np.float32);nlist.append('principalComponentScore152')
    global principalComponentScore153;principalComponentScore153=np.zeros(0);vlist.append(principalComponentScore153);dlist.append(np.float32);nlist.append('principalComponentScore153')
    global principalComponentScore154;principalComponentScore154=np.zeros(0);vlist.append(principalComponentScore154);dlist.append(np.float32);nlist.append('principalComponentScore154')
    global principalComponentScore155;principalComponentScore155=np.zeros(0);vlist.append(principalComponentScore155);dlist.append(np.float32);nlist.append('principalComponentScore155')
    global principalComponentScore156;principalComponentScore156=np.zeros(0);vlist.append(principalComponentScore156);dlist.append(np.float32);nlist.append('principalComponentScore156')
    global principalComponentScore157;principalComponentScore157=np.zeros(0);vlist.append(principalComponentScore157);dlist.append(np.float32);nlist.append('principalComponentScore157')
    global principalComponentScore158;principalComponentScore158=np.zeros(0);vlist.append(principalComponentScore158);dlist.append(np.float32);nlist.append('principalComponentScore158')
    global principalComponentScore159;principalComponentScore159=np.zeros(0);vlist.append(principalComponentScore159);dlist.append(np.float32);nlist.append('principalComponentScore159')
    global principalComponentScore160;principalComponentScore160=np.zeros(0);vlist.append(principalComponentScore160);dlist.append(np.float32);nlist.append('principalComponentScore160')
    global principalComponentScore161;principalComponentScore161=np.zeros(0);vlist.append(principalComponentScore161);dlist.append(np.float32);nlist.append('principalComponentScore161')
    global principalComponentScore162;principalComponentScore162=np.zeros(0);vlist.append(principalComponentScore162);dlist.append(np.float32);nlist.append('principalComponentScore162')
    global principalComponentScore163;principalComponentScore163=np.zeros(0);vlist.append(principalComponentScore163);dlist.append(np.float32);nlist.append('principalComponentScore163')
    global principalComponentScore164;principalComponentScore164=np.zeros(0);vlist.append(principalComponentScore164);dlist.append(np.float32);nlist.append('principalComponentScore164')
    global principalComponentScore165;principalComponentScore165=np.zeros(0);vlist.append(principalComponentScore165);dlist.append(np.float32);nlist.append('principalComponentScore165')
    global principalComponentScore166;principalComponentScore166=np.zeros(0);vlist.append(principalComponentScore166);dlist.append(np.float32);nlist.append('principalComponentScore166')
    global principalComponentScore167;principalComponentScore167=np.zeros(0);vlist.append(principalComponentScore167);dlist.append(np.float32);nlist.append('principalComponentScore167')
    global principalComponentScore168;principalComponentScore168=np.zeros(0);vlist.append(principalComponentScore168);dlist.append(np.float32);nlist.append('principalComponentScore168')
    global principalComponentScore169;principalComponentScore169=np.zeros(0);vlist.append(principalComponentScore169);dlist.append(np.float32);nlist.append('principalComponentScore169')
    global principalComponentScore170;principalComponentScore170=np.zeros(0);vlist.append(principalComponentScore170);dlist.append(np.float32);nlist.append('principalComponentScore170')
    global principalComponentScore171;principalComponentScore171=np.zeros(0);vlist.append(principalComponentScore171);dlist.append(np.float32);nlist.append('principalComponentScore171')
    global principalComponentScore172;principalComponentScore172=np.zeros(0);vlist.append(principalComponentScore172);dlist.append(np.float32);nlist.append('principalComponentScore172')
    global principalComponentScore173;principalComponentScore173=np.zeros(0);vlist.append(principalComponentScore173);dlist.append(np.float32);nlist.append('principalComponentScore173')
    global principalComponentScore174;principalComponentScore174=np.zeros(0);vlist.append(principalComponentScore174);dlist.append(np.float32);nlist.append('principalComponentScore174')
    global principalComponentScore175;principalComponentScore175=np.zeros(0);vlist.append(principalComponentScore175);dlist.append(np.float32);nlist.append('principalComponentScore175')
    global principalComponentScore176;principalComponentScore176=np.zeros(0);vlist.append(principalComponentScore176);dlist.append(np.float32);nlist.append('principalComponentScore176')
    global principalComponentScore177;principalComponentScore177=np.zeros(0);vlist.append(principalComponentScore177);dlist.append(np.float32);nlist.append('principalComponentScore177')
    global principalComponentScore178;principalComponentScore178=np.zeros(0);vlist.append(principalComponentScore178);dlist.append(np.float32);nlist.append('principalComponentScore178')
    global principalComponentScore179;principalComponentScore179=np.zeros(0);vlist.append(principalComponentScore179);dlist.append(np.float32);nlist.append('principalComponentScore179')
    global principalComponentScore180;principalComponentScore180=np.zeros(0);vlist.append(principalComponentScore180);dlist.append(np.float32);nlist.append('principalComponentScore180')
    global principalComponentScore181;principalComponentScore181=np.zeros(0);vlist.append(principalComponentScore181);dlist.append(np.float32);nlist.append('principalComponentScore181')
    global principalComponentScore182;principalComponentScore182=np.zeros(0);vlist.append(principalComponentScore182);dlist.append(np.float32);nlist.append('principalComponentScore182')
    global principalComponentScore183;principalComponentScore183=np.zeros(0);vlist.append(principalComponentScore183);dlist.append(np.float32);nlist.append('principalComponentScore183')
    global principalComponentScore184;principalComponentScore184=np.zeros(0);vlist.append(principalComponentScore184);dlist.append(np.float32);nlist.append('principalComponentScore184')
    global principalComponentScore185;principalComponentScore185=np.zeros(0);vlist.append(principalComponentScore185);dlist.append(np.float32);nlist.append('principalComponentScore185')
    global principalComponentScore186;principalComponentScore186=np.zeros(0);vlist.append(principalComponentScore186);dlist.append(np.float32);nlist.append('principalComponentScore186')
    global principalComponentScore187;principalComponentScore187=np.zeros(0);vlist.append(principalComponentScore187);dlist.append(np.float32);nlist.append('principalComponentScore187')
    global principalComponentScore188;principalComponentScore188=np.zeros(0);vlist.append(principalComponentScore188);dlist.append(np.float32);nlist.append('principalComponentScore188')
    global principalComponentScore189;principalComponentScore189=np.zeros(0);vlist.append(principalComponentScore189);dlist.append(np.float32);nlist.append('principalComponentScore189')
    global principalComponentScore190;principalComponentScore190=np.zeros(0);vlist.append(principalComponentScore190);dlist.append(np.float32);nlist.append('principalComponentScore190')
    global principalComponentScore191;principalComponentScore191=np.zeros(0);vlist.append(principalComponentScore191);dlist.append(np.float32);nlist.append('principalComponentScore191')
    global principalComponentScore192;principalComponentScore192=np.zeros(0);vlist.append(principalComponentScore192);dlist.append(np.float32);nlist.append('principalComponentScore192')
    global principalComponentScore193;principalComponentScore193=np.zeros(0);vlist.append(principalComponentScore193);dlist.append(np.float32);nlist.append('principalComponentScore193')
    global principalComponentScore194;principalComponentScore194=np.zeros(0);vlist.append(principalComponentScore194);dlist.append(np.float32);nlist.append('principalComponentScore194')
    global principalComponentScore195;principalComponentScore195=np.zeros(0);vlist.append(principalComponentScore195);dlist.append(np.float32);nlist.append('principalComponentScore195')
    global principalComponentScore196;principalComponentScore196=np.zeros(0);vlist.append(principalComponentScore196);dlist.append(np.float32);nlist.append('principalComponentScore196')
    global principalComponentScore197;principalComponentScore197=np.zeros(0);vlist.append(principalComponentScore197);dlist.append(np.float32);nlist.append('principalComponentScore197')
    global principalComponentScore198;principalComponentScore198=np.zeros(0);vlist.append(principalComponentScore198);dlist.append(np.float32);nlist.append('principalComponentScore198')
    global principalComponentScore199;principalComponentScore199=np.zeros(0);vlist.append(principalComponentScore199);dlist.append(np.float32);nlist.append('principalComponentScore199')
    global principalComponentScore200;principalComponentScore200=np.zeros(0);vlist.append(principalComponentScore200);dlist.append(np.float32);nlist.append('principalComponentScore200')
    global principalComponentScore201;principalComponentScore201=np.zeros(0);vlist.append(principalComponentScore201);dlist.append(np.float32);nlist.append('principalComponentScore201')
    global principalComponentScore202;principalComponentScore202=np.zeros(0);vlist.append(principalComponentScore202);dlist.append(np.float32);nlist.append('principalComponentScore202')
    global principalComponentScore203;principalComponentScore203=np.zeros(0);vlist.append(principalComponentScore203);dlist.append(np.float32);nlist.append('principalComponentScore203')
    global principalComponentScore204;principalComponentScore204=np.zeros(0);vlist.append(principalComponentScore204);dlist.append(np.float32);nlist.append('principalComponentScore204')
    global principalComponentScore205;principalComponentScore205=np.zeros(0);vlist.append(principalComponentScore205);dlist.append(np.float32);nlist.append('principalComponentScore205')
    global principalComponentScore206;principalComponentScore206=np.zeros(0);vlist.append(principalComponentScore206);dlist.append(np.float32);nlist.append('principalComponentScore206')
    global principalComponentScore207;principalComponentScore207=np.zeros(0);vlist.append(principalComponentScore207);dlist.append(np.float32);nlist.append('principalComponentScore207')
    global principalComponentScore208;principalComponentScore208=np.zeros(0);vlist.append(principalComponentScore208);dlist.append(np.float32);nlist.append('principalComponentScore208')
    global principalComponentScore209;principalComponentScore209=np.zeros(0);vlist.append(principalComponentScore209);dlist.append(np.float32);nlist.append('principalComponentScore209')
    global principalComponentScore210;principalComponentScore210=np.zeros(0);vlist.append(principalComponentScore210);dlist.append(np.float32);nlist.append('principalComponentScore210')
    global principalComponentScore211;principalComponentScore211=np.zeros(0);vlist.append(principalComponentScore211);dlist.append(np.float32);nlist.append('principalComponentScore211')
    global principalComponentScore212;principalComponentScore212=np.zeros(0);vlist.append(principalComponentScore212);dlist.append(np.float32);nlist.append('principalComponentScore212')
    global principalComponentScore213;principalComponentScore213=np.zeros(0);vlist.append(principalComponentScore213);dlist.append(np.float32);nlist.append('principalComponentScore213')
    global principalComponentScore214;principalComponentScore214=np.zeros(0);vlist.append(principalComponentScore214);dlist.append(np.float32);nlist.append('principalComponentScore214')
    global principalComponentScore215;principalComponentScore215=np.zeros(0);vlist.append(principalComponentScore215);dlist.append(np.float32);nlist.append('principalComponentScore215')
    global principalComponentScore216;principalComponentScore216=np.zeros(0);vlist.append(principalComponentScore216);dlist.append(np.float32);nlist.append('principalComponentScore216')
    global principalComponentScore217;principalComponentScore217=np.zeros(0);vlist.append(principalComponentScore217);dlist.append(np.float32);nlist.append('principalComponentScore217')
    global principalComponentScore218;principalComponentScore218=np.zeros(0);vlist.append(principalComponentScore218);dlist.append(np.float32);nlist.append('principalComponentScore218')
    global principalComponentScore219;principalComponentScore219=np.zeros(0);vlist.append(principalComponentScore219);dlist.append(np.float32);nlist.append('principalComponentScore219')
    global principalComponentScore220;principalComponentScore220=np.zeros(0);vlist.append(principalComponentScore220);dlist.append(np.float32);nlist.append('principalComponentScore220')
    global principalComponentScore221;principalComponentScore221=np.zeros(0);vlist.append(principalComponentScore221);dlist.append(np.float32);nlist.append('principalComponentScore221')
    global principalComponentScore222;principalComponentScore222=np.zeros(0);vlist.append(principalComponentScore222);dlist.append(np.float32);nlist.append('principalComponentScore222')
    global principalComponentScore223;principalComponentScore223=np.zeros(0);vlist.append(principalComponentScore223);dlist.append(np.float32);nlist.append('principalComponentScore223')
    global principalComponentScore224;principalComponentScore224=np.zeros(0);vlist.append(principalComponentScore224);dlist.append(np.float32);nlist.append('principalComponentScore224')
    global principalComponentScore225;principalComponentScore225=np.zeros(0);vlist.append(principalComponentScore225);dlist.append(np.float32);nlist.append('principalComponentScore225')
    global principalComponentScore226;principalComponentScore226=np.zeros(0);vlist.append(principalComponentScore226);dlist.append(np.float32);nlist.append('principalComponentScore226')
    global principalComponentScore227;principalComponentScore227=np.zeros(0);vlist.append(principalComponentScore227);dlist.append(np.float32);nlist.append('principalComponentScore227')
    global principalComponentScore228;principalComponentScore228=np.zeros(0);vlist.append(principalComponentScore228);dlist.append(np.float32);nlist.append('principalComponentScore228')
    global principalComponentScore229;principalComponentScore229=np.zeros(0);vlist.append(principalComponentScore229);dlist.append(np.float32);nlist.append('principalComponentScore229')
    global principalComponentScore230;principalComponentScore230=np.zeros(0);vlist.append(principalComponentScore230);dlist.append(np.float32);nlist.append('principalComponentScore230')
    global principalComponentScore231;principalComponentScore231=np.zeros(0);vlist.append(principalComponentScore231);dlist.append(np.float32);nlist.append('principalComponentScore231')
    global principalComponentScore232;principalComponentScore232=np.zeros(0);vlist.append(principalComponentScore232);dlist.append(np.float32);nlist.append('principalComponentScore232')
    global principalComponentScore233;principalComponentScore233=np.zeros(0);vlist.append(principalComponentScore233);dlist.append(np.float32);nlist.append('principalComponentScore233')
    global principalComponentScore234;principalComponentScore234=np.zeros(0);vlist.append(principalComponentScore234);dlist.append(np.float32);nlist.append('principalComponentScore234')
    global principalComponentScore235;principalComponentScore235=np.zeros(0);vlist.append(principalComponentScore235);dlist.append(np.float32);nlist.append('principalComponentScore235')
    global principalComponentScore236;principalComponentScore236=np.zeros(0);vlist.append(principalComponentScore236);dlist.append(np.float32);nlist.append('principalComponentScore236')
    global principalComponentScore237;principalComponentScore237=np.zeros(0);vlist.append(principalComponentScore237);dlist.append(np.float32);nlist.append('principalComponentScore237')
    global principalComponentScore238;principalComponentScore238=np.zeros(0);vlist.append(principalComponentScore238);dlist.append(np.float32);nlist.append('principalComponentScore238')
    global principalComponentScore239;principalComponentScore239=np.zeros(0);vlist.append(principalComponentScore239);dlist.append(np.float32);nlist.append('principalComponentScore239')
    global principalComponentScore240;principalComponentScore240=np.zeros(0);vlist.append(principalComponentScore240);dlist.append(np.float32);nlist.append('principalComponentScore240')
    global principalComponentScore241;principalComponentScore241=np.zeros(0);vlist.append(principalComponentScore241);dlist.append(np.float32);nlist.append('principalComponentScore241')
    global principalComponentScore242;principalComponentScore242=np.zeros(0);vlist.append(principalComponentScore242);dlist.append(np.float32);nlist.append('principalComponentScore242')
    global principalComponentScore243;principalComponentScore243=np.zeros(0);vlist.append(principalComponentScore243);dlist.append(np.float32);nlist.append('principalComponentScore243')
    global principalComponentScore244;principalComponentScore244=np.zeros(0);vlist.append(principalComponentScore244);dlist.append(np.float32);nlist.append('principalComponentScore244')
    global principalComponentScore245;principalComponentScore245=np.zeros(0);vlist.append(principalComponentScore245);dlist.append(np.float32);nlist.append('principalComponentScore245')
    global principalComponentScore246;principalComponentScore246=np.zeros(0);vlist.append(principalComponentScore246);dlist.append(np.float32);nlist.append('principalComponentScore246')
    global principalComponentScore247;principalComponentScore247=np.zeros(0);vlist.append(principalComponentScore247);dlist.append(np.float32);nlist.append('principalComponentScore247')
    global principalComponentScore248;principalComponentScore248=np.zeros(0);vlist.append(principalComponentScore248);dlist.append(np.float32);nlist.append('principalComponentScore248')
    global principalComponentScore249;principalComponentScore249=np.zeros(0);vlist.append(principalComponentScore249);dlist.append(np.float32);nlist.append('principalComponentScore249')
    global principalComponentScore250;principalComponentScore250=np.zeros(0);vlist.append(principalComponentScore250);dlist.append(np.float32);nlist.append('principalComponentScore250')
    global principalComponentScore251;principalComponentScore251=np.zeros(0);vlist.append(principalComponentScore251);dlist.append(np.float32);nlist.append('principalComponentScore251')
    global principalComponentScore252;principalComponentScore252=np.zeros(0);vlist.append(principalComponentScore252);dlist.append(np.float32);nlist.append('principalComponentScore252')
    global principalComponentScore253;principalComponentScore253=np.zeros(0);vlist.append(principalComponentScore253);dlist.append(np.float32);nlist.append('principalComponentScore253')
    global principalComponentScore254;principalComponentScore254=np.zeros(0);vlist.append(principalComponentScore254);dlist.append(np.float32);nlist.append('principalComponentScore254')
    global principalComponentScore255;principalComponentScore255=np.zeros(0);vlist.append(principalComponentScore255);dlist.append(np.float32);nlist.append('principalComponentScore255')
    global principalComponentScore256;principalComponentScore256=np.zeros(0);vlist.append(principalComponentScore256);dlist.append(np.float32);nlist.append('principalComponentScore256')
    global principalComponentScore257;principalComponentScore257=np.zeros(0);vlist.append(principalComponentScore257);dlist.append(np.float32);nlist.append('principalComponentScore257')
    global principalComponentScore258;principalComponentScore258=np.zeros(0);vlist.append(principalComponentScore258);dlist.append(np.float32);nlist.append('principalComponentScore258')
    global principalComponentScore259;principalComponentScore259=np.zeros(0);vlist.append(principalComponentScore259);dlist.append(np.float32);nlist.append('principalComponentScore259')
    global principalComponentScore260;principalComponentScore260=np.zeros(0);vlist.append(principalComponentScore260);dlist.append(np.float32);nlist.append('principalComponentScore260')
    global principalComponentScore261;principalComponentScore261=np.zeros(0);vlist.append(principalComponentScore261);dlist.append(np.float32);nlist.append('principalComponentScore261')
    global principalComponentScore262;principalComponentScore262=np.zeros(0);vlist.append(principalComponentScore262);dlist.append(np.float32);nlist.append('principalComponentScore262')
    global principalComponentScore263;principalComponentScore263=np.zeros(0);vlist.append(principalComponentScore263);dlist.append(np.float32);nlist.append('principalComponentScore263')
    global principalComponentScore264;principalComponentScore264=np.zeros(0);vlist.append(principalComponentScore264);dlist.append(np.float32);nlist.append('principalComponentScore264')
    global principalComponentScore265;principalComponentScore265=np.zeros(0);vlist.append(principalComponentScore265);dlist.append(np.float32);nlist.append('principalComponentScore265')
    global principalComponentScore266;principalComponentScore266=np.zeros(0);vlist.append(principalComponentScore266);dlist.append(np.float32);nlist.append('principalComponentScore266')
    global principalComponentScore267;principalComponentScore267=np.zeros(0);vlist.append(principalComponentScore267);dlist.append(np.float32);nlist.append('principalComponentScore267')
    global principalComponentScore268;principalComponentScore268=np.zeros(0);vlist.append(principalComponentScore268);dlist.append(np.float32);nlist.append('principalComponentScore268')
    global principalComponentScore269;principalComponentScore269=np.zeros(0);vlist.append(principalComponentScore269);dlist.append(np.float32);nlist.append('principalComponentScore269')
    global principalComponentScore270;principalComponentScore270=np.zeros(0);vlist.append(principalComponentScore270);dlist.append(np.float32);nlist.append('principalComponentScore270')
    global principalComponentScore271;principalComponentScore271=np.zeros(0);vlist.append(principalComponentScore271);dlist.append(np.float32);nlist.append('principalComponentScore271')
    global principalComponentScore272;principalComponentScore272=np.zeros(0);vlist.append(principalComponentScore272);dlist.append(np.float32);nlist.append('principalComponentScore272')
    global principalComponentScore273;principalComponentScore273=np.zeros(0);vlist.append(principalComponentScore273);dlist.append(np.float32);nlist.append('principalComponentScore273')
    global principalComponentScore274;principalComponentScore274=np.zeros(0);vlist.append(principalComponentScore274);dlist.append(np.float32);nlist.append('principalComponentScore274')
    global principalComponentScore275;principalComponentScore275=np.zeros(0);vlist.append(principalComponentScore275);dlist.append(np.float32);nlist.append('principalComponentScore275')
    global principalComponentScore276;principalComponentScore276=np.zeros(0);vlist.append(principalComponentScore276);dlist.append(np.float32);nlist.append('principalComponentScore276')
    global principalComponentScore277;principalComponentScore277=np.zeros(0);vlist.append(principalComponentScore277);dlist.append(np.float32);nlist.append('principalComponentScore277')
    global principalComponentScore278;principalComponentScore278=np.zeros(0);vlist.append(principalComponentScore278);dlist.append(np.float32);nlist.append('principalComponentScore278')
    global principalComponentScore279;principalComponentScore279=np.zeros(0);vlist.append(principalComponentScore279);dlist.append(np.float32);nlist.append('principalComponentScore279')
    global principalComponentScore280;principalComponentScore280=np.zeros(0);vlist.append(principalComponentScore280);dlist.append(np.float32);nlist.append('principalComponentScore280')
    global principalComponentScore281;principalComponentScore281=np.zeros(0);vlist.append(principalComponentScore281);dlist.append(np.float32);nlist.append('principalComponentScore281')
    global principalComponentScore282;principalComponentScore282=np.zeros(0);vlist.append(principalComponentScore282);dlist.append(np.float32);nlist.append('principalComponentScore282')
    global principalComponentScore283;principalComponentScore283=np.zeros(0);vlist.append(principalComponentScore283);dlist.append(np.float32);nlist.append('principalComponentScore283')
    global principalComponentScore284;principalComponentScore284=np.zeros(0);vlist.append(principalComponentScore284);dlist.append(np.float32);nlist.append('principalComponentScore284')
    global principalComponentScore285;principalComponentScore285=np.zeros(0);vlist.append(principalComponentScore285);dlist.append(np.float32);nlist.append('principalComponentScore285')
    global principalComponentScore286;principalComponentScore286=np.zeros(0);vlist.append(principalComponentScore286);dlist.append(np.float32);nlist.append('principalComponentScore286')
    global principalComponentScore287;principalComponentScore287=np.zeros(0);vlist.append(principalComponentScore287);dlist.append(np.float32);nlist.append('principalComponentScore287')
    global principalComponentScore288;principalComponentScore288=np.zeros(0);vlist.append(principalComponentScore288);dlist.append(np.float32);nlist.append('principalComponentScore288')
    global principalComponentScore289;principalComponentScore289=np.zeros(0);vlist.append(principalComponentScore289);dlist.append(np.float32);nlist.append('principalComponentScore289')
    global principalComponentScore290;principalComponentScore290=np.zeros(0);vlist.append(principalComponentScore290);dlist.append(np.float32);nlist.append('principalComponentScore290')
    global principalComponentScore291;principalComponentScore291=np.zeros(0);vlist.append(principalComponentScore291);dlist.append(np.float32);nlist.append('principalComponentScore291')
    global principalComponentScore292;principalComponentScore292=np.zeros(0);vlist.append(principalComponentScore292);dlist.append(np.float32);nlist.append('principalComponentScore292')
    global principalComponentScore293;principalComponentScore293=np.zeros(0);vlist.append(principalComponentScore293);dlist.append(np.float32);nlist.append('principalComponentScore293')
    global principalComponentScore294;principalComponentScore294=np.zeros(0);vlist.append(principalComponentScore294);dlist.append(np.float32);nlist.append('principalComponentScore294')
    global principalComponentScore295;principalComponentScore295=np.zeros(0);vlist.append(principalComponentScore295);dlist.append(np.float32);nlist.append('principalComponentScore295')
    global principalComponentScore296;principalComponentScore296=np.zeros(0);vlist.append(principalComponentScore296);dlist.append(np.float32);nlist.append('principalComponentScore296')
    global principalComponentScore297;principalComponentScore297=np.zeros(0);vlist.append(principalComponentScore297);dlist.append(np.float32);nlist.append('principalComponentScore297')
    global principalComponentScore298;principalComponentScore298=np.zeros(0);vlist.append(principalComponentScore298);dlist.append(np.float32);nlist.append('principalComponentScore298')
    global principalComponentScore299;principalComponentScore299=np.zeros(0);vlist.append(principalComponentScore299);dlist.append(np.float32);nlist.append('principalComponentScore299')
    global principalComponentScore300;principalComponentScore300=np.zeros(0);vlist.append(principalComponentScore300);dlist.append(np.float32);nlist.append('principalComponentScore300')
    return vlist,dlist,nlist

if __name__ == "__main__":
    main()

