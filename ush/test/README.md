## Prerequisite
   1. Clone and build obsForge  

      ```       
      git clone --recursive https://github.com/noaa-emc/obsForge -b feature/initial 
      cd ./obsForge/sorc/ioda
      git remote -v
      git remote set-url origin https://github.com/jcsda-internal/ioda.git
      git pull
      git checkout feature/bufr_in_parallel_emily
      ./build.sh
      ```

   2. Clone wxflow (no need to build)

      ```
      git clone https://github.com/NOAA-EMC/wxflow 
      ```

   Example: obsForge and wxflow builds on HERA
   - obsForge: /scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge      
   - wxflow:   /scratch1/NCEPDEV/da/Emily.Liu/EMC-wxflow/wxflow
   

## Elements should be in the working directory (e.g. run_satwnd)
- Required input files:
   - bufr2ioda_bufr_backend_satwnd_amv_goes.yaml
   - bufr2ioda_satwnd_amv_goes_mapping.yaml
   - bufr2ioda_satwnd_amv_goes.py
   - bufr2ioda_script_backend_satwnd_amv_goes.yaml
   - testinput/2021080100/gdas.t00z.satwnd.tm00.bufr_d

- Processing shell script:
   - process_bufr2ioda

## How to run the test shell script
- Get help page for usage

```
      ./process_bufrioda -h

      <obsforge_dir> : root directory of obsForge build
      <wxflow_dir>   : root directory of wxflow build
      <cycle>        : cycle time (e.g., 2021080100)
      <obstype>      : observation type to create (e.g., satwnd_amv_goes)
      <sensor>       : sensor (e.g., abi)
      <mode>         : mode of operation (three valid modes: bufr_backend, script_backend, and bufr2netcdf)
      <nproc>        : number of processors (must be a positive integer)
```

- Run with default input parameters 

```
      ./process_bufrioda
```

- Run with user-defined input parameters 

```
      obsforge_dir="/scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge"
      wxflow_dir="/scratch1/NCEPDEV/da/Emily.Liu/EMC-wxflow/wxflow"
      ./process_bufr2ioda ${obsforge_dir} ${wxflow_dir} 2021080100 satwnd_amv_goes abi bufr2netcdf 12
```

-  Run with user-defined mode 

```
     ./process_bufr2ioda "" "" "" "" "" bufr2netcdf" 4 
     ./process_bufr2ioda "" "" "" "" "" bufr_backend" 4 
     ./process_bufr2ioda "" "" "" "" "" script_backend" 4 

```
