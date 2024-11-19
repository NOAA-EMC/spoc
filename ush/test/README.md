## Test bufr_query mapping, Python converter script, and ioda configuration YAML in obsForge
This is a prototype for testing BUFR to IODA conversion and is still evolving.

## Prerequisite
- Clone and build obsForge  

   ```       
      git clone --recursive https://github.com/noaa-emc/obsForge -b feature/initial
   
      cd ./obsForge/sorc/ioda
   
      git remote -v
   
      git remote set-url origin https://github.com/jcsda-internal/ioda.git
   
      git pull
   
      git checkout feature/bufr_in_parallel_emily
   
      ./build.sh
   ```

- Clone wxflow (no need to build)

   ```
      git clone https://github.com/NOAA-EMC/wxflow 
   ```

- Example: obsForge and wxflow builds on HERA
   ```
      obsForge  /scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge
        
      wxflow    /scratch1/NCEPDEV/da/Emily.Liu/EMC-wxflow/wxflow
   ```

## Elements should be in the working directory from SPOC
- Required input files:
  
   - bufr2ioda_satwnd_amv_goes.py
     
   - bufr2ioda_satwnd_amv_goes_mapping.yaml
     
   - bufr2ioda_bufr_backend_satwnd_amv_goes.yaml
     
   - bufr2ioda_script_backend_satwnd_amv_goes.yaml
     
   - testinput/2021080100/gdas.t00z.satwnd.tm00.bufr_d (copied from the global dump)

- Processing shell script:
   - bufr2ioda.sh 

## How to run the test shell script
- Get the help page for usage

```
      ./bufrioda.sh -h

      <obsforge_dir> : root directory of obsForge build
      <wxflow_dir>   : root directory of wxflow build
      <cycle>        : cycle time (e.g., 2021080100)
      <obstype>      : observation type to create (e.g., satwnd_amv_goes)
      <sensor>       : sensor (e.g., abi)
      <mode>         : mode of operation (four valid modes: bufr_backend, script_backend, bufr2netcdf, script2netcdf)
      <nproc>        : number of processors (must be a positive integer)
```

- Run with default input parameters 

```
      ./bufrioda.sh
```

- Run with user-defined input parameters 

```
      obsforge_dir="/scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge"

      wxflow_dir="/scratch1/NCEPDEV/da/Emily.Liu/EMC-wxflow/wxflow"

      ./bufr2ioda.sh ${obsforge_dir} ${wxflow_dir} 2021080100 satwnd_amv_goes abi script_backend 4 
```

-  Run with user-defined mode and number of processes

```
     ./bufr2ioda.sh "" "" "" "" "" bufr2netcdf" 8 

     ./bufr2ioda.sh "" "" "" "" "" script2netcdf" 0 

     ./bufr2ioda.sh "" "" "" "" "" bufr_backend" 12 

```
