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


def compute_orbit_direction(latitude, timestamp):
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


def compute_scan_position(fovn, cut_spot=11, dfov=4.25):
    """
    Compute scan position from Field of View Number (FOVN) for AVHRR GAC data.
    This is the python version of the gsi Fortran code 
    in src/gsi/read_avhrr.f90 that converts 
    FOVN from the bufr file to scan position:

    !
    !       Get scan position (1 - 90) based on (409 - 2*cut_spot - 1) = 386 here,  GAC pixels
    !       avhrr gac scan has 409 positions. we drop tails: 1- 11 & 399- 409 [the two ends]
    !       here we linearly map pixels: 12- 398 to 1 - 90 scan positions
               if ( mod(hdr(10)-cut_spot,dfov) < half*dfov ) then
                  scan_pos = real(nint((hdr(10)-cut_spot)/dfov) + 1) 
               else
                  scan_pos = real(nint((hdr(10)-cut_spot)/dfov)) 
               endif

    """
    half_dfov = 0.5 * dfov
    offset = fovn - cut_spot

    scan_position = ma.masked_array(np.zeros_like(fovn, dtype=float), mask=fovn.mask)

    condition = ma.masked_less((offset % dfov), half_dfov).filled(False)

    scan_position[condition] = np.round(offset[condition] / dfov) + 1
    scan_position[~condition] = np.round(offset[~condition] / dfov)

    return ma.masked_where(fovn.mask, scan_position.astype(int))


def compute_sensor_view_angle(scan_position, max_angle=55.4):
    """
    Compute sensor view angle from scan position.
    """
    num_positions = 90
    center = (num_positions + 1) / 2
    angle_per_pos = max_angle / (center - 1)

    deviation = scan_position - center
    view_angle = deviation * angle_per_pos
    return ma.masked_where(scan_position.mask, view_angle)


def compute_sensor_azimuth_angle(sensor_view_angle, orbit_direction):
    """
    Approximate sensor azimuth angle based on orbit direction.
    """
    azimuth = ma.masked_array(np.zeros_like(sensor_view_angle), mask=sensor_view_angle.mask)

    ascending = orbit_direction == 1
    descending = orbit_direction == -1

    azimuth[ascending] = 90 + sensor_view_angle[ascending]
    azimuth[descending] = 270 - sensor_view_angle[descending]

    return azimuth


def compute_solar_azimuth_angle(latitude, longitude, timestamp):
    """
    Approximate solar azimuth angle using a basic astronomical model.
    For accurate calculations, consider using `pvlib` or `skyfield`.
    """
    solar_azimuth = ma.masked_array(np.zeros_like(latitude), mask=latitude.mask)

    for i in range(len(latitude)):
        if latitude.mask[i] or longitude.mask[i] or timestamp.mask[i]:
            continue

        dt = datetime.utcfromtimestamp(timestamp[i])
        hour_angle = ((dt.hour + dt.minute / 60 + dt.second / 3600) - 12) * 15
        day_of_year = dt.timetuple().tm_yday
        decl = 23.44 * np.cos(np.deg2rad((360 / 365.0) * (day_of_year + 10)))
        lat = latitude[i]

        cos_az = (
            np.sin(np.deg2rad(decl)) - np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(90))
        ) / (np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(90)))
        
        az = np.rad2deg(np.arccos(cos_az))
        solar_azimuth[i] = az

    return solar_azimuth


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

            # get the data to compute additional variables
            timestamp = container.get('timestamp', cat)
            latitude = container.get('latitude', cat)
            longitude = container.get('longitude', cat)
            if latitude.size != 0:
                self.log.debug(f'latitude min/max = {latitude.min()} {latitude.max()}')
            if longitude.size != 0:
                self.log.debug(f'longitude min/max = {longitude.min()} {longitude.max()}')       
            if timestamp.size != 0:
                self.log.debug(f'timestamp min/max = {timestamp.min()} {timestamp.max()}')                                                                           
            fovn = container.get("fieldOfViewNumber", cat)

            # compute the additional variables
            scan_position = compute_scan_position(fovn)
            sensor_view_angle = compute_sensor_view_angle(scan_position)
            orbit_direction = compute_orbit_direction(latitude, timestamp)
            sensor_azimuth_angle = compute_sensor_azimuth_angle(sensor_view_angle, orbit_direction)
            solar_azimuth_angle = compute_solar_azimuth_angle(latitude, longitude, timestamp)

            # add the new variables to the container
            paths = container.get_paths("sensorZenithAngle", cat)
            container.add("sensorScanPosition", scan_position, paths, cat)
            # print(f'>>> scan_position type = {scan_position.dtype}')
            container.add("sensorViewAngle", sensor_view_angle, paths, cat)
            container.add("sensorAzimuthAngle", sensor_azimuth_angle, paths, cat)
            container.add("solarAzimuthAngle", solar_azimuth_angle, paths, cat)

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
