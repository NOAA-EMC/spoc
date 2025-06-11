#!/usr/bin/env python3

import os
import sys
import numpy as np
import numpy.ma as ma
import netCDF4 as nc

import bufr
from bufr.obs_builder import ObsBuilder


def map_path(map_file_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, map_file_name)

# OceanBasin class provides a facility to add an OceanBasin
# metadata variable using lon and lat
# basic definition of ocean basins is read from an nc file,
# We search for the filename, depending on the system
# The path to the ocean basin nc file can be supplied
# in the implementation of the converter

# the main method is get_station_basin which returns the ocean basin
# for a list of station coordinates

class OceanBasin:
    def __init__(self, file_path):
        self.latitudes = None
        self.longitudes = None
        self.basin_array = None

        self.read_nc_file(file_path)

    def read_nc_file(self, file_path):
        try:
            with nc.Dataset(file_path, 'r') as nc_file:
                variable_name = 'open_ocean'
                if variable_name in nc_file.variables:
                    lat_dim = nc_file.dimensions['lat'].size
                    lon_dim = nc_file.dimensions['lon'].size
                    self.latitudes = nc_file.variables['lat'][:]
                    self.longitudes = nc_file.variables['lon'][:]

                    variable = nc_file.variables[variable_name]
                    # Read the variable data into a numpy array
                    variable_data = variable[:]
                    # Convert to 2D numpy array
                    self.basin_array = np.reshape(variable_data, (lat_dim, lon_dim))
        except FileNotFoundError:
            print(f"The file {file_path} does not exist.")
            sys.exit(1)
        except IOError as e:
            # Handle other I/O errors, such as permission errors
            print(f"An IOError occurred: {e}")
            sys.exit(1)

    # input: 2 vectors of station coordinates
    # output: a vector of station ocean basin values
    def get_station_basin(self, lat, lon):
        n = len(lon)

        lat0 = self.latitudes[0]
        dlat = self.latitudes[1] - self.latitudes[0]
        lon0 = self.longitudes[0]
        dlon = self.longitudes[1] - self.longitudes[0]

        # the data may be a masked array
        # ocean_basin = []
        # for i in range(n):
            # if not ma.is_masked(lat[i]):
                # i1 = round((lat[i] - lat0) / dlat)
                # i2 = round((lon[i] - lon0) / dlon)
                # ocean_basin.append(self.basin_array[i1][i2])
        # return np.array(ocean_basin, dtype=np.int32)


        # the data may be a masked array
        ocean_basin = ma.array([0]*lat.size, mask=lat.mask, dtype=np.int32)
        for i in range(n):
            if not ma.is_masked(lat[i]):
                i1 = round((lat[i] - lat0) / dlat)
                i2 = round((lon[i] - lon0) / dlon)
                ocean_basin[i] = self.basin_array[i1][i2]
        return ocean_basin



def clean_lat_lon(lat, lon):
    """
    return a mask for valid pairs of latitude and longitude.

    Parameters:
    lat : array-like or None
        Latitude values, expected in [-90, 90].
    lon : array-like or None
        Longitude values, all in [-180, 180] or all in [0, 360].

    Returns:
    mask : numpy.ndarray
        Boolean mask (same shape as inputs), True for valid lat/lon pairs.
        Use to filter other arrays (e.g., other_data[mask]).
    """
    # Handle undefined or None inputs
    if lat is None or lon is None:
        return np.array([], dtype=bool)

    # Convert inputs to NumPy arrays, preserving masked arrays
    try:
        lat = np.asarray(lat, dtype=float)
        lon = np.asarray(lon, dtype=float)
    except (ValueError, TypeError):
        return np.zeros_like(lat, dtype=bool)

    # Ensure arrays have the same shape
    # probably needs to throw an exception
    if lat.shape != lon.shape:
        return np.zeros_like(lat, dtype=bool)

    # Initialize mask (True for valid, False for invalid)
    mask = np.ones(lat.shape, dtype=bool)

    # Handle masked arrays
    if np.ma.isMaskedArray(lat):
        mask &= ~lat.mask
    if np.ma.isMaskedArray(lon):
        mask &= ~lon.mask

    # Validate latitude: must be in [-90, 90]
    mask &= (lat >= -90) & (lat <= 90) & ~np.isnan(lat)

    # Validate longitude: all in [-180, 180] or all in [0, 360]
    valid_lons = lon[mask]  # Longitudes where lat is valid
    if len(valid_lons) == 0:
        return mask

    # Check longitude range consistency
    in_neg180_180 = (valid_lons >= -180) & (valid_lons <= 180)
    in_0_360 = (valid_lons >= 0) & (valid_lons <= 360)

    # Determine if all valid longitudes are in one range
    all_neg180_180 = np.all(in_neg180_180)
    all_0_360 = np.all(in_0_360)

    # If neither range is fully consistent, use dominant range
    # if the fraction of the "other" range is too large,
    # probably needs to throw an error
    if not (all_neg180_180 or all_0_360):
        if np.sum(in_neg180_180) >= np.sum(in_0_360):
            mask &= (lon >= -180) & (lon <= 180)
        else:
            mask &= (lon >= 0) & (lon <= 360)
    else:
        if all_neg180_180:
            mask &= (lon >= -180) & (lon <= 180)
        else:
            mask &= (lon >= 0) & (lon <= 360)

    # Ensure no NaN in longitude
    mask &= ~np.isnan(lon)
    return mask


class MarineInsituObsBuilder(ObsBuilder):
    def __init__(self, mapping_path, log_name=os.path.basename(__file__), config=None):
        self.ocean_basin_file = config['ocean_basin'] if 'ocean_basin' in config else None
        super().__init__(mapping_path, log_name=log_name, config=config)

    def make_obs(self, comm, input_path):
        container = super().make_obs(comm, input_path)

        lat = container.get('latitude')
        lon = container.get('longitude')
        lat_lon_mask = clean_lat_lon(lat, lon)
        container.apply_mask(lat_lon_mask)

        if self.ocean_basin_file and os.path.exists(self.ocean_basin_file):
            self._add_ocean_basin(container, self.ocean_basin_file)
        else:
            self.log.warning(f"No ocean basin file provided, or can not be found")
        return container

    def _make_description(self):
        description = super()._make_description()

        if self.ocean_basin_file and os.path.exists(self.ocean_basin_file):
            description.add_variables([
                {'name': "MetaData/oceanBasin",
                 'source': 'oceanBasin',
                 'longName': "Ocean Basin",
                 'units': ""}])

        return description


    def _add_preqc_var(self, container, name):
        v = container.get(name)
        paths = container.get_paths(name)
        preqc_name = f"PreQC{name}"
        preqc = np.zeros_like(v)
        container.add(preqc_name, preqc, paths)

    def _add_error_var(self, container, name, error):
        v = container.get(name)
        paths = container.get_paths(name)
        # print(f"_add_error_var {name}")
        error_var_name = f"ObsError{name}"
        # print(f"_add_error_var {error_var_name}")
        error_var = np.full_like(v, error)
        container.add(error_var_name, error_var, paths)

    def _add_seq_num(self, container):
        lon = container.get("longitude")
        lat = container.get("latitude")
        paths = container.get_paths("latitude")
        combined = np.stack((lon, lat), axis=-1)
        unique_combined, seq_num = np.unique(combined, axis=0, return_inverse=True)
        v = np.ma.masked_array(seq_num, mask=lon.mask)
        container.add("SequenceNumber", v, paths)

    def _add_ocean_basin(self, container, nc_file_path):
        lon = container.get("longitude")
        lat = container.get("latitude")
        paths = container.get_paths("latitude")
        ocean = OceanBasin(nc_file_path)
        v = ocean.get_station_basin(lat, lon)
        container.add("oceanBasin", v, paths)
