## Test bufr_query mapping, Python converter script, and ioda configuration YAML in obsForge
This is a prototype for testing BUFR to NETCDF and other formats (IODA, ZARR, ...etc.) conversion and is still evolving.

## Prerequisite
- Clone and build obsForge  

   ```       
      git clone --recursive https://github.com/noaa-emc/obsForge
   
      cd ./obsForge/sorc/ioda
   
      git remote -v
   
      git remote set-url origin https://github.com/jcsda-internal/ioda.git
   
      git pull
   
      git checkout feature/bufr_in_parallel

      cd ../../
   
      ./build.sh
   ```

- Example: obsForge builds on HERA
  
   ```
      obsForge  /scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge
   ```

## Elements should be in the working directory from SPOC
Creat a working directory (e.g. work_dir)

- Required input files in ./work_dir:
  
   - bufr_satwnd_amv_abi.py (copied /spoc/dump/mapping)
     
   - bufr_satwnd_amv_abi_mapping.yaml (copied from /spoc/dump/mapping)
     
   - bufr_bufr4backend_satwnd_amv_abi.yaml (copied from /spoc/ush/test/config)
     
   - bufr_script4backend_satwnd_amv_abi.yaml (copied from /spoc/ush/test/config)
     
   - /testinput/2021080100/gdas.t00z.satwnd.tm00.bufr_d (copied from the global dump)

- Processing shell script in ./work_dir :
   - ./encodeBufr.sh (copied from /spoc/ush/test)

## How to run the test shell script
- Get the help page for usage

```
      encodeBufr.sh -h

      <obsforge_dir>      : root directory of obsForge build
      <cycle>             : cycle time (e.g., 2021080100)
      <bufrtype>          : BUFR dump type to process (e.g., satwnd, atms, sfcsno)
      <obstype>           : observation type to create (e.g., satwnd_amv_abi, atms, sfcsno)
      <sensor>            : sensor (e.g., abi, atms); for non-satellite dta, sensor is usually obstype (e.g., sfcsno)
      <split_by_category> : split the data into multiple files based on category (false or true)
      <mode>              : mode of operation (e.g., bufr4backend, script4backend, bufr2netcdf, script2netcdf)
      <nproc>             : number of processors (positive integer to run with MPI, or zero for serial execution)
```

- Run with default input parameters 

```
      encodeBufr.sh
```

- Run with user-defined input parameters 

```
      obsforge_dir="/scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge"

      encodeBufr.sh ${obsforge_dir} 2021080100 satwnd satwnd_amv_abi abi true script4backend 4 

      encodeBufr.sh ${obsforge_dir} 2021080100 sfcsno sfcsno sfcsno false script4backend 4 

      encodeBufr.sh ${obsforge_dir} 2021080100 atms atms atms true script4backend 4 
```

-  Run with user-defined mode and number of processes

```
     encodeBufr.sh "" "" "" "" "" "" bufr2netcdf" 8 

     encodeBufr.sh "" "" "" "" "" "" script2netcdf" 0 

     encodeBufr.sh "" "" "" "" "" "" bufr4backend" 12 

     encodeBufr.sh "" "" "" "" "" "" script4backend" 4
```
