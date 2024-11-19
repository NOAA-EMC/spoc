#!/bin/bash

# ====================
# Get input arguments
# ====================
#obsforge_dir="${1:-/scratch1/NCEPDEV/da/Emily.Liu/EMC-obsForge/obsForge}"
obsforge_dir="${1:-/scratch1/NCEPDEV/da/Emily.Liu/obsforge}"
wxflow_dir="${2:-/scratch1/NCEPDEV/da/Emily.Liu/EMC-wxflow/wxflow}"
cycle="${3:-2021080100}"
obstype="${4:-satwnd_amv_goes}"
sensor="${5:-abi}"
mode="${6:-script_backend}"
nproc="${7:-4}"

# ==========================
# Function to display usage
# ==========================
usage() {
    echo "Usage: $0 <obsforge_dir> <wxflow> <cycle> <obstype> <sensor> <mode> <nproc>"
    echo "  <obsforge_dir> : root directory of obsForge build"
    echo "  <wxflow_dir>   : root directory of wxflow build"
    echo "  <cycle>        : cycle time (e.g., 2021080100)"
    echo "  <obstype>      : observation type to create (e.g., satwnd_amv_goes)"
    echo "  <sensor>       : sensor (e.g., abi)"
    echo "  <mode>         : mode of operation (e.g., bufr_backend, script_backend, bufr2netcdf, script2netcdf)"
    echo "  <nproc>        : number of processors (must be a positive integer to run with MPI or zero to run without MPI)"
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
if [[ -z "$mode" || -z "$nproc" || -z "$sensor" || -z "$obstype" || -z "$cycle" || -z "${obsforge_dir}" || -z "${wxflow_dir}" ]]; then
    echo "Error: Missing arguments."
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
echo "root directory of wxflow: $wxflow_dir"
echo "cycle: $cycle"
echo "obstype: $obstype"
echo "sensor: $sensor"
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

# =====================
# Set WXFLOW log level 
# =====================
export LOG_LEVEL=DEBUG

# ===============================
# Load obsForge required modules 
# ===============================
module use ${obsforge_dir}/modulefiles
module load obsforge/hera.intel
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
export PYTHONPATH="${PYTHONPATH}:${wxflow_dir}/src"

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
in_dir="${work_dir}/testinput/"
mkdir -p -m770 ${out_dir}

# ===============
# Set file paths
# ===============
mapping_file="${work_dir}/bufr2ioda_${obstype}_mapping.yaml"
input_file="${in_dir}/gdas.t${h2}z.satwnd.tm00.bufr_d"
output_file="${out_dir}/gdas.t${h2}z.satwnd_${sensor}_{splits/satId}.tm00.nc"
ioda_config_yaml="${work_dir}/bufr2ioda_${mode}_${obstype}.yaml"

# =============================
# Run ioda bufr/script backend
# =============================
if [[ "$mode" == "bufr_backend" || "$mode" == "script_backend" ]]; then
   if [[ "$nproc" == "0" ]]; then
      echo Run time_IodaIO.x without MPI ...
      ${obsforge_dir}/build/bin/time_IodaIO.x ${ioda_config_yaml}
   else
      echo Run time_IodaIO.x with MPI ${nproc} ...
      srun -n $nproc  --mem 96G --time 00:30:00 ${obsforge_dir}/build/bin/time_IodaIO.x ${ioda_config_yaml}
   fi
# ================
# Run bufr2netcdf  
# ================
elif [[ "$mode" == "bufr2netcdf" ]]; then
   if [[ "$nproc" == "0" ]]; then
      echo Run bufr2netcdf without MPI ...
      ${obsforge_dir}/build/bin/bufr2netcdf.x "$input_file" "${mapping_file}" "$output_file"
   else
      echo Run bufr2netcdf with MPI ${nproc} ...
      srun -n "$nproc" --time 00:30:00 --mem 96G ${obsforge_dir}/build/bin/bufr2netcdf.x "$input_file" "${mapping_file}" "$output_file" 
   fi
# ==================
# Run script2netcdf  
# ==================
elif [[ "$mode" == "script2netcdf" ]]; then
   if [[ "$nproc" == "0" ]]; then
      echo Run script2netcdf without MPI ...
      python bufr2ioda_${obstype}.py -m "$mapping_file" -o "$output_file" -i "$input_file"
   else
      echo Run script2netcdf with MPI ${nproc} ...
      srun -n "$nproc" --time 00:30:00 --mem 96G python bufr2ioda_${obstype}.py -m "$mapping_file" -o "$output_file" -i "$input_file"
   fi
else
   echo Incorrect running mode ${mode} ... Valid modes are: bufr_backend, script_back, bufr2netcdf, or script2netcdf
fi


