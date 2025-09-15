#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma
import sys
import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path
import pvlib
import wxflow
from bufr.transforms import compute_solar_angles


MAPPING_PATH = map_path('bufr_avhrr.yaml')


def determine_orbit_direction(latitude, timestamp):
    """
    Determine ascending or descending pass for each point.
    
    :param latitude: List or numpy array of latitude (degrees).
    :param timestamp: List or numpy array of UNIX timestamp (seconds).
    
    :return: List of 'ascending' or 'descending' strings corresponding to each interval.
    """
    directions = np.zeros_like(latitude)
    for i in range(1, len(latitude)):
        if latitude[i] > latitude[i-1]:
            directions[i] = 1
        elif latitude[i] < latitude[i-1]:
            directions[i] = -1
        else:
            directions[i] = 0
    return directions


def get_sensor_scan_position(fov_array):
    # Normalize FOV scan position to [-1, 1]
    scan_position = np.zeros_like(fov_array)
    scan_position = scan_position.astype(np.float32)
    if fov_array.size == 0:
        return scan_position

    min_fov = fov_array.min()
    max_fov = fov_array.max()
    scan_position = 2 * ((fov_array - min_fov) / (max_fov - min_fov)) - 1
    print(f'xxxxx    scan_position type = {scan_position.dtype}')
    scan_position = scan_position.astype(np.float32)
    print(f'yyyyy  scan_position type = {scan_position.dtype}')
    return scan_position


def compute_avhrr_sensor_azimuth(timestamp, latitude, fov_array):
    """
    Compute sensor azimuth angle for AVHRR data based on FOV and orbit direction.

    Parameters:
        timestamp   : 1D array of UNIX timestamp (same length as latitude)
        latitude    : 1D array of latitude in degrees
        fov_array    : 1D or masked array of field-of-view numbers (integers)

    Returns:
        azimuth      : Masked array of sensor azimuth angles in degrees [0, 360)
    """
    # Ensure masked array
    fov_array = ma.masked_array(fov_array)

    # Prepare output
    azimuth = ma.masked_all(fov_array.shape, dtype=np.float32)
    if fov_array.size == 0:
        return azimuth

    scan_position = get_sensor_scan_position(fov_array)

    # Get orbit direction: 1 = descending, -1 = ascending
    orbit_direction = determine_orbit_direction(latitude, timestamp)
    orbit_direction = np.array(orbit_direction)  # Ensure it's an array for masking

    # Descending pass (north to south)
    descending_mask = (orbit_direction == 1) & (~fov_array.mask)
    azimuth[descending_mask] = 180 + 90 * scan_position[descending_mask]

    # Ascending pass (south to north)
    ascending_mask = (orbit_direction == -1) & (~fov_array.mask)
    azimuth[ascending_mask] = 180 - 90 * scan_position[ascending_mask]

    # Wrap azimuth to [0, 360)
    azimuth = azimuth % 360

    return azimuth


def compute_sensor_scan_and_view_angle(fov_array, total_fovs=2048, max_view_angle_deg=55.0):
    """
    Compute sensorScanPosition and sensorViewAngle for AVHRR-style instruments.

    Parameters:
        fov_array : np.ma.MaskedArray
            Field-of-view number array (0 to total_fovs-1)
        total_fovs : int
            Total number of fields-of-view (default 2048)
        max_view_angle_deg : float
            Max scan angle from nadir (default 55.0° for AVHRR)

    Returns:
        Tuple of np.ma.MaskedArray:
            - sensorScanPosition (unitless, -1.0 to +1.0)
            - sensorViewAngle (degrees)
    """
    fov_array = ma.masked_array(fov_array)

    # Normalize scan position: -1 to +1
    scan_pos = 2 * (fov_array / (total_fovs - 1)) - 1

    # View angle in degrees (approximate)
    view_angle = scan_pos * max_view_angle_deg

    return scan_pos, view_angle


class BufrGsrasrObsBuilder(ObsBuilder):

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'Container list (original): {container.list()}')
        self.log.debug(f'All_sub_categories =  {container.all_sub_categories()}')
        self.log.debug(f'Category map = {container.get_category_map()}')

        for cat in container.all_sub_categories():
            self.log.debug(f'category: {cat}')

            temp = container.get("brightnessTemperature", cat)

            name = "brightnessTemperature"
            v = container.get(name, cat)
            paths = container.get_paths(name, cat)

            preqc_name = f"PreQC{name}"
            preqc = np.zeros_like(v)
            container.add(preqc_name, preqc, paths, cat)

            error_var_name = f"ObsError{name}"
            error = 0.1
            error_var = np.full_like(v, error)
            container.add(error_var_name, error_var, paths, cat)

            # compute sensorAzimuthAngle
            timestamp = container.get('timestamp', cat)
            latitude = container.get('latitude', cat)
            longitude = container.get('longitude', cat)
            if latitude.size != 0:
                self.log.debug(f'latitude min/max = {latitude.min()} {latitude.max()}')
            if longitude.size != 0:
                self.log.debug(f'longitude min/max = {longitude.min()} {longitude.max()}')       
            if timestamp.size != 0:
                self.log.debug(f'timestamp min/max = {timestamp.min()} {timestamp.max()}')                                                                           

            # solar azimuth angle
            solar_zenith_angle, solar_azimuth_angle = compute_solar_angles(latitude, longitude, timestamp)
            paths = container.get_paths("sensorZenithAngle", cat)
            azimuth_angle_name = "solarAzimuthAngle"
            container.add(azimuth_angle_name, solar_azimuth_angle, paths, cat)

            # sensor azimuth angle
            fovn = container.get("fieldOfViewNumber", cat)
            sensor_azimuth_angle = compute_avhrr_sensor_azimuth(timestamp, latitude, fovn)
            azimuth_angle_name = "sensorAzimuthAngle"
            container.add(azimuth_angle_name, sensor_azimuth_angle, paths, cat)

            # sensor scan position
            scan_position = get_sensor_scan_position(fovn)
            scan_position_name = "sensorScanPosition"
            container.add(scan_position_name, scan_position, paths, cat)
            print(f'>>> scan_position type = {scan_position.dtype}')

            # sensor view angle
            max_view_angle = 55.4       # Max view angle for AVHRR (approx)
            view_angle = scan_position * max_view_angle
            view_angle_name = "sensorViewAngle"
            container.add(view_angle_name, view_angle, paths, cat)

            satId = container.get('satelliteId', cat)

            if not np.any(satId):
                self.log.warning(f'Warning: empty category {cat[0]} in input file')
                continue

        # Final container state
        self.log.debug(f'Container list (updated): {container.list()}')

        return container

    def _make_description(self):
        description = super()._make_description()
        return description


add_main_functions(BufrGsrasrObsBuilder)
