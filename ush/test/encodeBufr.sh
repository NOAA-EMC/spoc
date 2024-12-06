#!/bin/bash

# ====================
# Get input arguments
# ====================
obsforge_dir="${1:-/scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge}"
cycle="${2:-2021080100}"
bufrtype="${3:-satwnd}"
obstype="${4:-satwnd_amv_goes}"
sensor="${5:-abi}"
split_by_category="${6:-true}"
mode="${7:-script_backend}"
nproc="${8:-4}"

# ==========================
# Function to display usage
# ==========================
usage() {
    echo "Usage: $0 <obsforge_dir> <cycle> <obstype> <sensor> <mode> <nproc>"
    echo "  <obsforge_dir>      : root directory of obsForge build"
    echo "  <cycle>             : cycle time (e.g., 2021080100)"
    echo "  <bufrtype>          : BUFR dump type to process (e.g., satwnd, atms, cris, sfcsno)"
    echo "  <obstype>           : observation type to create (e.g., satwnd_amv_goes, atms, cris, sfcsno )"
    echo "  <sensor>            : sensor (e.g., abi, atms); for non-satellite dta, sensor is usually obstype (e.g., sfcsno)"
    echo "  <split_by_category> : split the output file into multiple files based on category (false or true)"
    echo "  <mode>              : mode of operation (e.g., bufr_backend, script_backend, bufr2netcdf, script2netcdf)"
    echo "  <nproc>             : number of processors (positive integer to run with MPI, or zero for serial execution)"
    exit 1
}

# =============================
# Check for --help or -h flags
# =============================
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
fi

# =============================================
# Check if all required arguments are provided
# =============================================
if [[ -z "$mode" || -z "$nproc" || -z "$sensor" || -z "${split_by_category}" || -z "$obstype" || -z "$cycle" || -z "${obsforge_dir}" ]]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# ===============
# Validate cycle 
# ===============
if ! [[ "$cycle" =~ ^[0-9]{10}$ ]]; then
    echo "Error: Invalid cycle format. Expected YYYYMMDDHH."
    usage
fi

# ==============
# Validate mode
# ==============
if [[ "$mode" != "bufr_backend" && "$mode" != "script_backend" && "$mode" != "bufr2netcdf"  && "$mode" != "script2netcdf" ]]; then
    echo "Error: Invalid mode '$mode'. Expected 'bufr_backend' or 'script_backend' or 'bufr2netcdf' or 'script2netcdf'."
    usage
fi

# =============================================
# Validate nproc is a positive integer or zero
# =============================================
if ! [[ "$nproc" =~ ^[0-9]+$ && "$nproc" -ge 0 ]]; then
	echo "Error: nproc must be a positive integer or zero."
    usage
fi

# =========================
# Check if all checks pass
# =========================
echo "root directory of obsForge: $obsforge_dir"
echo "cycle: $cycle"
echo "bufrtype: $bufrtype"
echo "obstype: $obstype"
echo "sensor: $sensor"
echo "split_by_category: ${split_by_category}"
echo "mode: $mode"
echo "number of processors: $nproc"

# ======================
# Enable debugging mode
# ======================
set -x

# =========================
# Set unlimited stack size
# =========================
ulimit -s unlimited
ulimit -a

# ====================================
# Set OOPS run time output parameters 
# ====================================
export OOPS_TRACE=1
export OOPS_DEBUG=1

# =======================================================
# Set WXFLOW log level
# Valid Log_LEVEL: DEBUG, INFO, WARNING, ERROR, CRITICAL
# =======================================================
export LOG_LEVEL=DEBUG

# ===============================
# Load obsForge required modules 
# ===============================
module use ${obsforge_dir}/modulefiles
module load obsforge/hera.intel || { echo "Error loading obsforge module"; exit 1; }
module list

# ==============================
# Set bufr-query python library
# ==============================
export LD_LIBRARY_PATH="${obsforge_dir}/build/lib:${LD_LIBRARY_PATH}"
export PYTHONPATH="${PYTHONPATH}:${obsforge_dir}/build/lib/python3.10/site-packages"

# ========================
# Set ioda python library
# =========================
export PYTHONPATH="${PYTHONPATH}:${obsforge_dir}/build/lib/python3.10"

# ============
# Set wxfloww 
# ============
export PYTHONPATH="${PYTHONPATH}:${obsforge_dir}/sorc/wxflow/src"

# ===========================
# Configure SLUM environment
# ===========================
export SLURM_ACCOUNT=da-cpu
export SALLOC_ACCOUNT=$SLURM_ACCOUNT
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
export SLURM_QOS=debug

# ===================================
# Extract year, month, day, and hour
# ===================================
y4="${cycle:0:4}"
m2="${cycle:4:2}"
d2="${cycle:6:2}"
h2="${cycle:8:2}"

# ====================
# Set directory paths
# ====================
work_dir="$PWD"
out_dir="${work_dir}/testoutput/$cycle/${mode}"
in_dir="${work_dir}/testinput/$cycle"
mkdir -p -m770 ${out_dir} || { echo "Error creating output directory: ${out_dir}"; exit 1; }

# ===============
# Set file paths
# ===============
ioda_config_yaml="${work_dir}/bufr2ioda_${mode}_${obstype}.yaml"
mapping_file="${work_dir}/bufr2ioda_${obstype}_mapping.yaml"
input_file="${in_dir}/gdas.t${h2}z.${bufrtype}.tm00.bufr_d"
if [[ "${split_by_category}" = "true" ]]; then
   output_file="${out_dir}/gdas.t${h2}z.${bufrtype}.${sensor}_{splits/satId}.tm00.nc"
else
   output_file="${out_dir}/gdas.t${h2}z.${bufrtype}.tm00.nc"
fi

if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file not found: $input_file"
    exit 1
fi

if [[ ! -f "$mapping_file" ]]; then
    echo "Error: mapping file not found: $mapping"
    exit 1
fi

# =============================
# Run ioda bufr/script backend
# =============================
if [[ "$mode" == "bufr_backend" || "$mode" == "script_backend" ]]; then
   if [[ ! -f "$ioda_config_yaml" ]]; then
      echo "Error: ioda configuration file not found: $ioda_config_yaml"
      exit 1
   fi
   if [[ "$nproc" == "0" ]]; then
      echo Run time_IodaIO.x without MPI ...
      ${obsforge_dir}/build/bin/time_IodaIO.x ${ioda_config_yaml} || { echo "Error: time_IodaIO.x failed"; exit 1; }
   else
      echo Run time_IodaIO.x with MPI ${nproc} ...
      srun -n $nproc --mem 96G --time 00:30:00 ${obsforge_dir}/build/bin/time_IodaIO.x ${ioda_config_yaml} || { echo "Error: MPI time_IodaIO.x failed"; exit 1; }
   fi
# ================
# Run bufr2netcdf  
# ================
elif [[ "$mode" == "bufr2netcdf" ]]; then
   if [[ "$nproc" == "0" ]]; then
      echo Run bufr2netcdf without MPI ...
      ${obsforge_dir}/build/bin/bufr2netcdf.x "$input_file" "${mapping_file}" "$output_file" || { echo "Error: bufr2netcdf.x failed"; exit 1; }
   else
      echo Run bufr2netcdf with MPI ${nproc} ...
      srun -n "$nproc" --mem 96G --time 00:30:00 ${obsforge_dir}/build/bin/bufr2netcdf.x "$input_file" "${mapping_file}" "$output_file" || { echo "Error: MPI bufr2netcdf.x failed"; exit 1; }
   fi
# ==================
# Run script2netcdf  
# ==================
elif [[ "$mode" == "script2netcdf" ]]; then
   if [[ "$nproc" == "0" ]]; then
      echo Run script2netcdf without MPI ...
      python bufr2ioda_${obstype}.py "$input_file" "$mapping_file" "$output_file" || { echo "Error: Python script2netcdf failed"; python bufr2ioda_${obstype}.py --help; exit 1; }
   else
      echo Run script2netcdf with MPI ${nproc} ...
      srun -n "$nproc" --mem 96G --time 00:30:00 python bufr2ioda_${obstype}.py "$input_file" "$mapping_file" "$output_file" || { echo "Error: MPI Python script2netcdf failed"; python bufr2ioda_${obstype}.py --help; exit 1; }
   fi
else
   echo Incorrect running mode ${mode} ... Valid modes are: bufr_backend, script_back, bufr2netcdf, or script2netcdf
fi


