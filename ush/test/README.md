## Test bufr_query mapping, Python converter script, and ioda configuration YAML in obsForge
This is a prototype for testing BUFR to IODA conversion and is still evolving.

## Prerequisite
- Clone and build obsForge  

   ```       
      git clone --recursive https://github.com/noaa-emc/obsForge
   
      cd ./obsForge/sorc/ioda
   
      git remote -v
   
      git remote set-url origin https://github.com/jcsda-internal/ioda.git
   
      git pull
   
      git checkout feature/bufr_in_parallel

      cd ../spoc

      git checkout feature/dump_satwind_goes

      cd ../../
   
      ./build.sh
   ```

- Example: obsForge builds on HERA
  
   ```
      obsForge  /scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge
   ```

## Elements should be in the working directory from SPOC
- Required input files:
  
   - bufr2ioda_satwnd_amv_goes.py
     
   - bufr2ioda_satwnd_amv_goes_mapping.yaml
     
   - bufr2ioda_bufr_backend_satwnd_amv_goes.yaml
     
   - bufr2ioda_script_backend_satwnd_amv_goes.yaml
     
   - testinput/2021080100/gdas.t00z.satwnd.tm00.bufr_d (copied from the global dump)

- Processing shell script:
   - encodeBufr.sh 

## How to run the test shell script
- Get the help page for usage

```
      bufrioda.sh -h

      <obsforge_dir>      : root directory of obsForge build
      <cycle>             : cycle time (e.g., 2021080100)
      <bufrtype>          : BUFR dump type to process (e.g., satwnd, atms, sfcsno)
      <obstype>           : observation type to create (e.g., satwnd_amv_goes, atms, sfcsno)
      <sensor>            : sensor (e.g., abi, atms); for non-satellite dta, sensor is usually obstype (e.g., sfcsno)
      <split_by_category> : split the data into multiple files based on category (false or true)
      <mode>              : mode of operation (e.g., bufr_backend, script_backend, bufr2netcdf, script2netcdf)
      <nproc>             : number of processors (positive integer to run with MPI, or zero for serial execution)
```

- Run with default input parameters 

```
      encodeBufr.sh
```

- Run with user-defined input parameters 

```
      obsforge_dir="/scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge"

      encodeBufr.sh ${obsforge_dir} 2021080100 satwnd satwnd_amv_goes abi true script_backend 4 

      encodeBufr.sh ${obsforge_dir} 2021080100 sfcsno sfcsno sfcsno false script_backend 4 

      encodeBufr.sh ${obsforge_dir} 2021080100 atms atms atms true script_backend 4 
```

-  Run with user-defined mode and number of processes

```
     encodeBufr.sh "" "" "" "" "" "" bufr2netcdf" 8 

     encodeBufr.sh "" "" "" "" "" "" script2netcdf" 0 

     encodeBufr.sh "" "" "" "" "" "" bufr_backend" 12 

     encodeBufr.sh "" "" "" "" "" "" script_backend" 4
```
