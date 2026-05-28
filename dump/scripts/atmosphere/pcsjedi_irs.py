#!/usr/bin/env python3
"""
Convert MTG-IRS NetCDF data to IODA format.

Reads Meteosat Third Generation - Infrared Sounder dwell files,
selects high-quality pixels, extracts radiometric and principal component data,
and encodes to IODA format for data assimilation.
"""

import os
import sys
import bufr
import argparse
import time as tm
import numpy as np
import concurrent.futures
from netCDF4 import Dataset
from itertools import repeat
from bufr.encoders import netcdf
from datetime import datetime, timezone


def main():
    """Main entry point for MTG-IRS to IODA conversion."""

    # Timing variables for performance tracking
    time_read = 0
    time_concat = 0
    time_encode = 0

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Convert MTG-IRS NetCDF files to IODA format"
    )
    parser.add_argument("--input", dest="inpdir", type=str, required=True,
                        help="input files directory")
    parser.add_argument("--output", dest="iodout", type=str, required=True,
                        help="name of ioda output file")
    parser.add_argument("--yaml", dest="yaml", type=str, required=True,
                        help="path/filename of yaml file")
    parser.add_argument("--grid", dest="grid", type=int, default=160,
                        help="input inst grid size (default: 160)")
    parser.add_argument("--thin", dest="thin", type=int, default=32,
                        help="input thin grid size (default: 32)")
    parser.add_argument("--pick", dest="pick", type=int, default=6,
                        help="input pick grid size (default: 6)")
    parser.add_argument("--npcs", dest="npcs", type=int, default=150,
                        help="number of principal components (default: 150)")
    parser.add_argument("--nwrk", dest="nwrk", type=int, default=20,
                        help="number of workers (default: 20)")
    args = parser.parse_args()

    inpdir = args.inpdir
    iodout = args.iodout
    yaml_file = args.yaml
    grid = args.grid
    thin = args.thin
    pick = args.pick
    npcs = args.npcs
    nwrk = args.nwrk

    # Validate input directory
    if not os.path.isdir(inpdir):
        sys.stderr.write(f"Error: input directory '{inpdir}' does not exist\n")
        sys.exit(1)

    # Build list of input files
    files = checkfiles(inpdir)
    nfils = len(files)
    print(f'Number of input files: {nfils}')

    if nfils == 0:
        sys.stderr.write(f"Error: no NetCDF files found in '{inpdir}'\n")
        sys.exit(1)

    # Set up file ranges for parallel processing
    chunk = (nfils + nwrk - 1) // nwrk  # Ceiling division
    fbeg = [idx * chunk for idx in range(nwrk)]
    fend = [min(idx * chunk + chunk, nfils) for idx in range(nwrk)]

    # Create variable metadata lists
    vlist, dlist, nlist, nvar = varlist(npcs)

    # Parallel loop through MTG-IRS dwell files
    print(f"Processing {nfils} files with {nwrk} workers...")
    rtime = tm.time()

    with concurrent.futures.ProcessPoolExecutor(max_workers=nwrk) as executor:
        results = list(executor.map(
            readloop,
            fbeg, fend,
            repeat(inpdir),
            repeat(grid),
            repeat(thin),
            repeat(pick),
            repeat(npcs),
            repeat(nwrk),
            repeat(nvar)
        ))

    time_read = tm.time() - rtime
    rtime = tm.time()

    # Concatenate arrays from parallel segments
    total_reports = 0
    for npic, data_segment in enumerate(results):
        kpic = np.count_nonzero(data_segment[0, :])
        if kpic == 0:
            continue

        vpic = 0
        for aray in range(len(vlist)):
            vlist[aray] = np.concatenate((vlist[aray], data_segment[vpic, 0:kpic]))
            vpic += 1
        total_reports += kpic

    # Ensure correct data types
    vpic = 0
    for aray in range(len(vlist)):
        vlist[aray] = vlist[aray].astype(dlist[vpic])
        vpic += 1

    time_concat = tm.time() - rtime
    rtime = tm.time()

    # Create data container and encode to IODA
    container = bufr.DataContainer()
    description = bufr.encoders.Description(yaml_file)

    # Add metadata and principal component variables
    npic = 0
    for aray in range(len(vlist)):
        container.add(nlist[npic], vlist[aray], ['*'])
        npic += 1

    # Add channel data (wavenumber and channel numbers)
    sensorCentralWavenumber, sensorChannelNumber = channel_data(files[0])
    container.add('sensorCentralWavenumber', sensorCentralWavenumber, ['CHANNEL'])
    container.add('sensorChannelNumber', sensorChannelNumber, ['CHANNEL'])

    # Encode container to IODA NetCDF file
    netcdf.Encoder(description).encode(container, iodout)

    time_encode = tm.time() - rtime

    # Print timing summary
    print()
    print("=" * 50)
    print("PROCESSING SUMMARY")
    print("=" * 50)
    print(f"File reading & selection:   {time_read:8.2f} sec")
    print(f"Data concatenation:         {time_concat:8.2f} sec")
    print(f"IODA encoding:              {time_encode:8.2f} sec")
    print(f"Total processing time:      {time_read + time_concat + time_encode:8.2f} sec")
    print(f"Total reports processed:    {total_reports:8d}")
    print("=" * 50)
    print()


def readloop(fbeg, fend, inpdir, grid, thin, pick, npcs, nwork, nvar):
    """
    Read and process a range of MTG-IRS dwell files in parallel.

    Args:
        fbeg: Starting file index
        fend: Ending file index
        inpdir: Input directory path
        grid: Grid size
        thin: Thinning factor
        pick: Pick box size
        npcs: Number of principal components
        nwork: Number of workers
        nvar: Number of metadata variables

    Returns:
        data: numpy array containing extracted observations
    """

    # Get file list and thinning parameters
    files = checkfiles(inpdir)
    kmax, ibox, jbox = thinparm(grid, thin, pick)
    ibox = np.array(ibox)
    jbox = np.array(jbox)

    # Set up storage arrays
    chunk = (len(files) + nwork - 1) // nwork
    data = np.zeros((nvar + 2 * npcs, kmax * chunk), dtype=np.float32)

    mpic = 0
    print(f"Worker processing files {fbeg} to {fend-1} from {inpdir}")

    for fidx in range(fbeg, fend):
        if fidx >= len(files):
            break

        try:
            with Dataset(files[fidx]) as irs:
                data_segment = _process_dwell_file(
                    irs, ibox, jbox, kmax, grid, thin, npcs, nvar, mpic, data
                )
                if data_segment is not None:
                    processed_count, updated_mpic = data_segment
                    mpic = updated_mpic
        except (OSError, KeyError) as e:
            print(f"Warning: skipping file {files[fidx]}: {e}")
            continue

    return data


def _process_dwell_file(irs, ibox, jbox, kmax, grid, thin, npcs, nvar, mpic, data):
    """
    Extract and process data from a single MTG-IRS dwell file.

    Args:
        irs: netCDF Dataset object
        ibox, jbox: Inner box indices
        kmax: Maximum number of pixels per dwell
        grid: Grid size
        thin: Thinning factor
        npcs: Number of principal components
        nvar: Number of metadata variables
        mpic: Current position in output array
        data: Output data array (modified in place)

    Returns:
        Tuple of (processed_count, updated_mpic) or None if no valid data
    """

    # Build composite quality score
    overall_quality = np.zeros([grid, grid], dtype='int32')
    
    try:
        lw_quality = {}
        for k in irs['data/lwir/quality_band'].variables.keys():
            if 'warning' in k and 'number' not in k:
                lw_quality[k] = irs['data/lwir/quality_band/' + k][:]
        
        mw_quality = {}
        for k in irs['data/mwir/quality_band'].variables.keys():
            if 'warning' in k and 'number' not in k:
                mw_quality[k] = irs['data/mwir/quality_band/' + k][:]
        
        for k in mw_quality.keys():
            overall_quality += (mw_quality[k].astype('int32') + 
                               lw_quality[k].astype('int32'))
    except KeyError:
        print("Warning: quality band data not found")

    # Extract key variables
    try:
        loca = irs['data']
        mwva = irs['data/mwir/compressed']
        lwva = irs['data/lwir/compressed']

        pcsc = lwva.variables["global_pc_scores"][:, :, :]
        lons = loca.variables["longitude"][:, :]
        lats = loca.variables["latitude"][:, :]
        satz = loca.variables["satellite_zenith_angle"][:, :]
    except KeyError as e:
        print(f"Warning: required variable missing: {e}")
        return None

    # Apply QC filters
    fil = 2147483647  # Max int32
    qc1 = np.abs(pcsc[:, :, :]) < fil
    qc2_lon = np.abs(lons) <= 180
    qc2_lat = np.abs(lats) <= 90
    qc2_zen = np.abs(satz) <= 60
    qc2 = qc2_lon & qc2_lat & qc2_zen

    pcsc = np.where(qc1, pcsc[:, :, :], 0)
    pcsc[:, :, 0] = np.where(qc2, pcsc[:, :, 0], 0)

    # Get satellite altitude
    try:
        sat_alt = irs['state/platform/platform_altitude'][0]
    except KeyError:
        sat_alt = 0.0

    # Select hottest spot in each inner box
    imax_list = []
    jmax_list = []
    
    for a in range(0, grid, thin):
        for b in range(0, grid, thin):
            pc1 = np.abs(pcsc[a + ibox, b + jbox, 0])
            pcm = np.max(pc1)
            if pcm > 0:
                inx = np.where(pc1 == pcm)[0]
                imax_list.append(ibox[inx[0]])
                jmax_list.append(jbox[inx[0]])

    # Skip if no valid pixels found
    if len(imax_list) == 0:
        return None

    imax = np.array(imax_list, dtype=int)
    jmax = np.array(jmax_list, dtype=int)
    kpic = len(imax)
    lpic = mpic
    mpic_new = mpic + kpic

    # Extract and store data
    try:
        data[0, lpic:mpic_new] = loca.variables["time"][:]
        data[1, lpic:mpic_new] = loca.variables["dwell_number"][:]
        data[2, lpic:mpic_new] = loca.variables["dwell_type"][:]
        data[3, lpic:mpic_new] = np.array(lwva.variables["detector_sample_quality"])[imax, jmax]
        data[4, lpic:mpic_new] = np.array(mwva.variables["detector_sample_quality"])[imax, jmax]
        data[5, lpic:mpic_new] = np.array(lwva.variables["global_pcrs_quality"])[imax, jmax]
        data[6, lpic:mpic_new] = np.array(mwva.variables["global_pcrs_quality"])[imax, jmax]
        data[7, lpic:mpic_new] = np.array(lwva.variables["spatial_sample_quality"])[imax, jmax]
        data[8, lpic:mpic_new] = np.array(mwva.variables["spatial_sample_quality"])[imax, jmax]
        data[9, lpic:mpic_new] = overall_quality[imax, jmax]
        data[10, lpic:mpic_new] = assign_WMO_ID(irs.platform)
        data[11, lpic:mpic_new] = scan_pos(data[1, lpic:mpic_new], imax, jmax, kpic, grid)
        data[12, lpic:mpic_new] = np.array(loca.variables["cloud_fraction"])[imax, jmax]
        data[13, lpic:mpic_new] = np.array(loca.variables["cloud_signal"])[imax, jmax]
        data[14, lpic:mpic_new] = np.array(loca.variables["dust_warning"])[imax, jmax]
        data[15, lpic:mpic_new] = np.array(loca.variables["latitude"])[imax, jmax]
        data[16, lpic:mpic_new] = np.array(loca.variables["longitude"])[imax, jmax]
        data[17, lpic:mpic_new] = np.array(loca.variables["satellite_azimuth_angle"])[imax, jmax]
        data[18, lpic:mpic_new] = np.array(loca.variables["satellite_zenith_angle"])[imax, jmax]
        data[19, lpic:mpic_new] = np.array(loca.variables["solar_azimuth_angle"])[imax, jmax]
        data[20, lpic:mpic_new] = np.array(loca.variables["solar_zenith_angle"])[imax, jmax]
        data[21, lpic:mpic_new] = compute_scan_angle(np.array(loca.variables["satellite_zenith_angle"])[imax, jmax], sat_alt)  
        data[22, lpic:mpic_new] = np.array(lwva.variables["global_pcr_scores"])[imax, jmax]
        data[23, lpic:mpic_new] = np.array(mwva.variables["global_pcr_scores"])[imax, jmax]

        # Extract principal component scores
        lw_scores = np.array(lwva.variables["global_pc_scores"])[imax, jmax, :npcs]
        mw_scores = np.array(mwva.variables["global_pc_scores"])[imax, jmax, :npcs]

        data[nvar:nvar + npcs, lpic:mpic_new] = lw_scores.T
        data[nvar + npcs:nvar + 2 * npcs, lpic:mpic_new] = mw_scores.T

    except (KeyError, IndexError) as e:
        print(f"Warning: error extracting data: {e}")
        return None

    return (kpic, mpic_new)


def assign_WMO_ID(platform):
    """
    Convert platform name to WMO satellite identifier.

    Args:
        platform: Platform name (e.g., 'MTS1', 'MTS2')

    Returns:
        WMO satellite ID (int)
    """
    wmo_mapping = {
        'MTS1': 72,
        'MTS2': 75,
    }
    
    wmo_id = wmo_mapping.get(platform)
    if wmo_id is None:
        print(f"Warning: unknown satellite platform: {platform}")
        wmo_id = -999  # Missing value indicator
    
    return wmo_id


def compute_scan_angle(sensor_zenith, sensor_altitude, qc_flag=None):
    """
    Compute satellite scan angle from altitude and zenith angle.

    Uses the relationship: γ = arcsin((R / (R + h)) * sin(θ))
    where R is Earth radius, h is altitude, θ is zenith angle.

    Args:
        sensor_altitude: Sensor altitude in m
        sensor_zenith: Satellite zenith angle in degrees
        qc_flag: Optional QC flag array (default: all good)

    Returns:
        Scan angle in degrees
    """
    earth_mean_radius_km = 6378.1370  # WGS84

    d2r = np.pi / 180.0
    r2d = 180.0 / np.pi

    ratio = np.empty_like(sensor_zenith)

    # Initialize QC flag if not provided
    if qc_flag is None:
        qc_flag = np.zeros_like(sensor_zenith)

    # Compute scan angle for good data
    good = qc_flag == 0
    if np.sum(good) > 0:
        ratio[good] = earth_mean_radius_km / (earth_mean_radius_km + sensor_altitude / 1000.0)

    scanang = np.arcsin(ratio * np.sin(np.abs(sensor_zenith) * d2r)) * r2d

    return scanang


def scan_pos(dwell_num, imax, jmax, kpic, grid):
    """
    Compute linearized scan position from dwell number and pixel indices.

    Args:
        dwell_num: Array of dwell numbers
        imax: Array of row indices
        jmax: Array of column indices
        kpic: Number of pixels
        grid: Grid size

    Returns:
        Array of scan positions
    """
    scanpos = np.zeros(kpic, dtype=np.int32)
    for k in range(kpic):
        scanpos[k] = (
            np.int64(dwell_num[k]) * np.int64(25600)
            + np.int64(imax[k]) * np.int64(grid)
            + np.int64(jmax[k])
        )
    return scanpos


def checkfiles(inpdir):
    """
    Get list of readable NetCDF files in directory.

    Args:
        inpdir: Input directory path

    Returns:
        List of full file paths to NetCDF files

    Raises:
        SystemExit if directory is invalid or no files found
    """
    try:
        file_list = os.listdir(inpdir)
    except OSError as e:
        sys.stderr.write(f"Error reading input directory '{inpdir}': {e}\n")
        sys.exit(1)

    netcdf_files = [f for f in file_list if f.lower().endswith('.nc')]
    if not netcdf_files:
        sys.stderr.write(f"No NetCDF files found in input directory '{inpdir}'\n")
        sys.exit(1)

    file_keep = []
    for filename in netcdf_files:
        file_path = os.path.join(inpdir, filename)
        file_keep.append(file_path)

    if not file_keep:
        sys.stderr.write(f"No readable NetCDF files found in '{inpdir}'\n")
        sys.exit(1)

    return file_keep


def thinparm(grid, outer, inner):
    """
    Calculate thinning parameters and box indices.

    Args:
        grid: Total grid size
        outer: Outer thinning box size
        inner: Inner selection box size

    Returns:
        Tuple of (kmax, ibox, jbox) where kmax is number of boxes
        and ibox/jbox are lists of indices for inner box selection

    Raises:
        SystemExit if box sizes are invalid
    """
    if grid % outer != 0:
        print(f'Error: outer box size {outer} does not divide grid size {grid}')
        sys.exit(1)

    if inner < 1 or inner > outer:
        print(f'Error: inner and outer box sizes must be >= 1')
        sys.exit(1)

    kmax = (grid // outer) ** 2

    beg = outer // 2 - inner // 2
    end = beg + inner

    if end - beg != inner:
        print(f'Error: invalid inner box size {inner} for outer box {outer}')
        sys.exit(1)

    ibox = []
    jbox = []
    for i in range(beg, end):
        for j in range(beg, end):
            ibox.append(i)
            jbox.append(j)

    return kmax, ibox, jbox


def varlist(npcs=150):
    """
    Define metadata variables and principal component variable lists.

    Args:
        npcs: Number of principal components (default: 150)

    Returns:
        Tuple of (vlist, dlist, nlist, nvar) where:
        - vlist: list of initialized numpy arrays
        - dlist: list of data types
        - nlist: list of variable names
        - nvar: number of metadata variables
    """
    # Metadata variable configuration: (name, dtype)
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

    vlist = []
    dlist = []
    nlist = []

    # Add metadata variables
    for name, dtype in meta_config:
        nlist.append(name)
        dlist.append(dtype)
        vlist.append(np.zeros(0, dtype=dtype))

    # Add principal component score variables
    for i in range(1, 2 * npcs + 1):
        name = f'principalComponentScore{i}'
        nlist.append(name)
        dlist.append(np.float32)
        vlist.append(np.zeros(0, dtype=np.float32))

    nvar = len(meta_config)
    return vlist, dlist, nlist, nvar


def channel_data(file):
    """
    Extract channel wavenumbers and channel numbers from MTG-IRS file.

    Args:
        file: Path to MTG-IRS NetCDF file

    Returns:
        Tuple of (wavenumber_array, channel_number_array)
    """
    try:
        with Dataset(file) as irs:
            wn_lw = np.asarray(irs['data/lwir/wavenumber'][:]).astype('float32')
            wn_mw = np.asarray(irs['data/mwir/wavenumber'][:]).astype('float32')
            waves = np.concatenate([wn_lw, wn_mw])
            chans = np.arange(1, len(waves) + 1, dtype='int32')
    except (OSError, KeyError) as e:
        print(f"Erroi: error reading channel data: {e}")
        sys.exit(1)

    return waves, chans


if __name__ == "__main__":
    main()

