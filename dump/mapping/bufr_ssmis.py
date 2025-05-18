#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma
from multiprocessing import Pool, cpu_count

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path, nprocs_per_task
from bufr.transforms import compute_solar_angles


MAPPING_PATH = map_path('bufr_ssmis.yaml')


class BufrSsmisObsBuilder(ObsBuilder):
    """
    Class for building observations from SSMIS BUFR data.

    This class extends `ObsBuilder` to include specific logic for processing
    SSMIS data such as solar angles and satellite ascending/descending orbits.

    :param mapping_path: Path to the mapping file.
    :type mapping_path: str
    """

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):

        # Get container from mapping file first
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'container list (original): {container.list()}')
        self.log.debug(f'all_sub_categories =  {container.all_sub_categories()}')
        self.log.debug(f'category map =  {container.get_category_map()}')

        # Add new/derived data into container
        for cat in container.all_sub_categories():

            self.log.debug(f'category = {cat}')

            satId = container.get('satelliteId', cat)
            if not np.any(satId):
                self.log.warning(f'category {cat[0]} does not exist in input file')

            self._add_solar_angles(container, cat)
            self._add_satellite_ascend_descent_orbit(container, cat)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')
        self.log.debug(f'all_sub_categories {container.all_sub_categories()}')

        return container

    def _make_description(self):
        description = super()._make_description()
        self._add_new_variable_descriptions(description)

        return description

    def _add_new_variable_descriptions(self, description):
        description.add_variables([
            {
                'name': 'MetaData/satelliteAscendingFlag',
                'source': 'satelliteAscendingFlag',
                'units': '1',
                'longName': 'Satellite Ascending/Descending Orbit Flag (Ascend:1; Descend:-1)',
            },
            {
                'name': 'MetaData/solarZenithAngle',
                'source': 'solarZenithAngle',
                'units': 'degree',
                'longName': 'Solar Zenith Angle',
            },
            {
                'name': 'MetaData/solarAzimuthAngle',
                'source': 'solarAzimuthAngle',
                'units': 'degree',
                'longName': 'Solar Azimuth Angle',
            }])

    def _add_satellite_ascend_descent_orbit(self, container, category):
        """
        Determine satellite orbit type (ascending or descending) based on latitude changes.

        :param container: Observation data container.
        :type container: Container
        :param category: Data category to process.
        :type category: str
        """

        satId = container.get('satelliteId', category)

        if not satId.size:
            paths = container.get_paths('fieldOfViewNumber', category)
            dummy = container.get('fieldOfViewNumber', category)
            container.add('satelliteAscendingFlag', dummy, paths, category)
            return

        # Get data from container
        # ephemeris data - latitude values in order of time
        first_lat = container.get('latitude1', category)
        self.log.debug(f'first_lat min/max = {first_lat.min()} {first_lat.max()}')
        second_lat = container.get('latitude2', category)
        self.log.debug(f'second_lat min/max = {second_lat.min()} {second_lat.max()}')
        fovn = container.get('fieldOfViewNumber', category)
        self.log.debug(f'fovn min/max = {fovn.min()} {fovn.max()}')

        # Determine ascending/descending mode
        # Compare latitude between the first and second records
        orbit = np.where(second_lat > first_lat, 1, -1).astype(np.int32)

        self.log.debug(f'orbit min/max = {orbit.min()} {orbit.max()}')

        paths = container.get_paths('fieldOfViewNumber', category)
        self.log.debug(f'paths = {paths}')
        container.add('satelliteAscendingFlag', orbit, paths, category)

    def _compute_solar_angles_parallel(self, latitudes, longitudes, unix_times, nprocs=None):
        """
        Compute solar zenith and azimuth angles in parallel using multiprocessing.

        :param latitudes: Array of latitudes in degrees.
        :type latitudes: numpy.ndarray
        :param longitudes: Array of longitudes in degrees.
        :type longitudes: numpy.ndarray
        :param unix_times: Array of Unix timestamps (seconds since 1970-01-01T00:00:00Z).
        :type unix_times: numpy.ndarray
        :param nprocs: Number of processes to use. Defaults to the number of CPU cores.
        :type nprocs: int, optional

        :return: Two arrays: zenith angles and azimuth angles.
        :rtype: tuple(numpy.ndarray, numpy.ndarray)
        """

        assert len(latitudes) == len(longitudes) == len(unix_times), "Input arrays must be the same length"

        if nprocs is None:
            nprocs = nprocs_per_task()

        args_list = list(zip(latitudes, longitudes, unix_times))

        self.log.debug(f'Using {nprocs} processes to compute solar angles.')

        with Pool(nprocs) as pool:
            results = pool.starmap(compute_solar_angles, args_list)

        self.log.debug(f'Using {nprocs} processes to compute solar angles --- done')
        zenith_angles, azimuth_angles = zip(*results)
        zenith_angles = np.array(zenith_angles)
        azimuth_angles = np.array(azimuth_angles)

        return zenith_angles, azimuth_angles

    def _add_solar_angles(self, container, category):
        """
        Compute and add solar zenith and azimuth angles to the observation container.

        :param container: Observation data container.
        :type container: Container
        :param category: Data category to process.
        :type category: str
        """

        satId = container.get('satelliteId', category)
        if not satId.size:
            paths = container.get_paths('latitude', category)
            dummy = container.get('latitude', category)
            container.add('solarZenithAngle', dummy, paths, category)
            container.add('solarAzimuthAngle', dummy, paths, category)
            return

        # Prepare input arrays
        unix_times = container.get('timestamp', category)
        latitudes = container.get('latitude', category)
        longitudes = container.get('longitude', category)
        self.log.debug(f'latitudes min/max = {latitudes.min()} {latitudes.max()}')
        self.log.debug(f'longitudes min/max = {longitudes.min()} {longitudes.max()}')
        self.log.debug(f'unix_times min/max = {unix_times.min()} {unix_times.max()}')

        # Calculate solar angles
        zenith_angles, azimuth_angles = self._compute_solar_angles_parallel(latitudes, longitudes, unix_times)

        self.log.debug(f'zenith_angles min/max = {zenith_angles.min()} {zenith_angles.max()}')
        self.log.debug(f'azimuth_angles min/max = {azimuth_angles.min()} {azimuth_angles.max()}')

        # Add solar angles to contained
        paths = container.get_paths('latitude', category)
        self.log.debug(f'paths = {paths}')
        container.add('solarZenithAngle', zenith_angles, paths, category)
        container.add('solarAzimuthAngle', azimuth_angles, paths, category)


# Add main functions create_obs_file or create_obs_group
add_main_functions(BufrSsmisObsBuilder)
