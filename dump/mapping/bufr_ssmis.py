#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma
import pytz
from pysolar.solar import get_altitude, get_azimuth
from datetime import datetime, timezone
from multiprocessing import Pool, cpu_count

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions


def map_path(map_file_name):
    """
    Generate the full path to the mapping file.

    :param map_file_name: Name of the mapping file.
    :type map_file_name: str
    :return: Full path to the mapping file.
    :rtype: str
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, map_file_name)


def get_default_nprocs(default=8):
    """
    Determine the number of processes to use from environment variables or fallback to a default.

    Environment Variables:
        - `SLURM_CPUS_PER_TASK`: Auto-detected from Slurm.
        - `CPUS_PER_TASK`: Input provided by the user.

    Examples:
        1. Run with 1 Slurm task and reserve 12 CPU cores for this one task.
           This spawns 12 worker processes via `multiprocessing.Pool`.
           (SLURM_CPUS_PER_TASK = 12)
           Usage:
               srun -n 1 --cpus_per_tasks=12 python bufr_ssmis.py

        2. Run without Slurm, specifying the number of worker processes.
           Usage:
               export CPUS_PER_TASK=12
               python bufr_ssmis.py

    :param default: Fallback default if environment variables are not set. Defaults to 8.
    :type default: int
    :return: Number of processes to use.
    :rtype: int
    """
    env_nprocs = os.getenv("CPUS_PER_TASK")
    slurm_nprocs = os.environ.get("SLURM_CPUS_PER_TASK")

    print(f"CPUS_PER_TASK={env_nprocs}, SLURM_CPUS_PER_TASK={slurm_nprocs}")

    try:
        return int(env_nprocs or slurm_nprocs or default)
    except ValueError:
        print(f"[WARN] Invalid environment variable value. Falling back to default: {default}")
        return default


def compute_solar_angles(lat, lon, unix_time):
    """
    Compute solar zenith and azimuth angles for a geographical point.

    :param lat: Latitude in degrees.
    :type lat: float
    :param lon: Longitude in degrees.
    :type lon: float
    :param unix_time: Unix timestamp (seconds since 1970-01-01T00:00:00Z).
    :type unix_time: int
    :return: Tuple containing zenith angle and azimuth angle (in degrees).
    :rtype: tuple(float, float)
    """
    dt = datetime.fromtimestamp(int(unix_time), tz=timezone.utc)
    altitude = get_altitude(lat, lon, dt)
    azimuth = get_azimuth(lat, lon, dt)
    zenith = 90.0 - altitude

    return zenith, azimuth


MAPPING_PATH = map_path('bufr_ssmis.yaml')


class BufrSsmisObsBuilder(ObsBuilder):
    """
    Builder class for processing BUFR SSMIS observations.

    This class extends the `ObsBuilder` and provides methods to construct
    observation data, add derived variables, and manage mapping configurations.
    """

    def __init__(self):
        """
        Initialize the BufrSsmisObsBuilder with the mapping path and log name.
        """
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        """
        Create observations and add derived data to the observation container.

        :param comm: Communication object for distributed processing.
        :type comm: obj
        :param input_path: Path to the input file.
        :type input_path: str
        :return: Updated observation container.
        :rtype: obj
        """
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'container list (original): {container.list()}')
        self.log.debug(f'all_sub_categories = {container.all_sub_categories()}')
        self.log.debug(f'category map = {container.get_category_map()}')

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

    def _add_satellite_ascend_descent_orbit(self, container, category):
        """
        Determine satellite orbit type (ascending or descending) based on latitude changes.

        If the latitude increases, the orbit is ascending (flag = 1).
        If the latitude decreases, the orbit is descending (flag = -1).

        :param container: Observation data container.
        :type container: obj
        :param category: Observation category.
        :type category: obj
        """
        satId = container.get('satelliteId', category)

        if not satId.size:
            paths = container.get_paths('fieldOfViewNumber', category)
            dummy = container.get('fieldOfViewNumber', category)
            container.add('satelliteAscendingFlag', dummy, paths, category)
            return

        # Get data from container
        first_lat = container.get('latitude1', category)
        second_lat = container.get('latitude2', category)
        fovn = container.get('fieldOfViewNumber', category)

        # Determine ascending/descending mode
        orbit = np.where(second_lat > first_lat, 1, -1).astype(np.int32)

        paths = container.get_paths('fieldOfViewNumber', category)
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
        :param nprocs: Number of processes to use. If None, defaults to the number of CPU cores.
        :type nprocs: int, optional
        :return: Two arrays containing zenith angles and azimuth angles, respectively.
        :rtype: tuple(numpy.ndarray, numpy.ndarray)
        """
        assert len(latitudes) == len(longitudes) == len(unix_times), "Input arrays must be the same length"

        if nprocs is None:
            nprocs = get_default_nprocs()

        args_list = list(zip(latitudes, longitudes, unix_times))

        with Pool(nprocs) as pool:
            results = pool.starmap(compute_solar_angles, args_list)

        zenith_angles, azimuth_angles = zip(*results)
        return np.array(zenith_angles), np.array(azimuth_angles)


add_main_functions(BufrSsmisObsBuilder)
