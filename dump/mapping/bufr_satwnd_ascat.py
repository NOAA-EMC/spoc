#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path
from bufr.transforms import compute_wind_components


MAPPING_PATH = map_path('bufr_satwnd_ascat.yaml')


class BufrAscatObsBuilder(ObsBuilder):
    """
    A builder class to generate satellite wind observations from ASCAT BUFR input.

    Inherits from :class:`bufr.obs_builder.ObsBuilder` and uses a mapping file to
    extract wind speed and direction, compute wind vector components, and attach
    observation metadata.
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

                paths = container.get_paths('windSpeedAt10M', cat)
                dummy = container.get('windSpeedAt10M', cat)
                container.add('windEastward', dummy, paths, cat)
                container.add('windNorthward', dummy, paths, cat)

                paths = container.get_paths('satelliteId', cat)
                dummy = container.get('satelliteId', cat)
                container.add('obstype_windEastward', dummy, paths, cat)
                container.add('obstype_windNorthward', dummy, paths, cat)
                continue

            wdir = container.get('windDirectionAt10M', cat)
            wspd = container.get('windSpeedAt10M', cat)
            self.log.debug(f'wdir min/max = {wdir.min()} {wdir.max()}')
            self.log.debug(f'wspd min/max = {wspd.min()} {wspd.max()}')

            uob, vob = compute_wind_components(wspd, wdir)
            self.log.debug(f'uob min/max = {uob.min()} {uob.max()}')
            self.log.debug(f'vob min/max = {vob.min()} {vob.max()}')

            paths = container.get_paths('windSpeedAt10M', cat)
            container.add('windEastward', uob, paths, cat)
            container.add('windNorthward', vob, paths, cat)

            obstype = self._get_obs_type(container, cat)
            paths = container.get_paths('satelliteId', cat)
            container.add('obstype_windEastward', obstype, paths, cat)
            container.add('obstype_windNorthward', obstype, paths, cat)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')
        self.log.debug(f'all_sub_categories {container.all_sub_categories()}')

        return container

    def _get_obs_type(self, container, category):
        """
        Retrieve observation type values for wind components.

        :param container: Observation container from BUFR input.
        :type container: bufr.DataContainer
        :param category: Subcategory identifier (e.g., sensor/platform group).
        :type category: tuple

        :returns: Array filled with observation type (290) for the given category.
        :rtype: numpy.ndarray
        """

        satId = container.get('satelliteId', category)

        if not satId.size:
            paths = container.get_paths('satelliteId', cat)
            dummy = container.get('satelliteId', cat)
            container.add('obstype_windEastward', dummy, paths, cat)
            container.add('obstype_windNorthward', dummy, paths, cat)
            return

        obstype = np.full_like(satId, 290)

        return obstype

    def _make_description(self):
        description = super()._make_description()
        self._add_new_variable_descriptions(description)

        return description

    def _add_new_variable_descriptions(self, description):
        description.add_variables([
            {
                'name': 'ObsValue/windEastward',
                'source': 'windEastward',
                'units': 'm s-1',
                'longName': '10-meter U-Wind Component',
            },
            {
                'name': 'ObsValue/windNorthward',
                'source': 'windNorthward',
                'units': 'm s-1',
                'longName': '10-meter V-Wind Component',
            },
            {
                'name': 'ObsType/windEastward',
                'source': 'obstype_windEastward',
                'units': '1',
                'longName': 'Observation Type for Wind Components',
            },
            {
                'name': 'ObsType/windNorthward',
                'source': 'obstype_windNorthward',
                'units': '1',
                'longName': 'Observation Type for Wind Components',
            }])


# Add main functions create_obs_file or create_obs_group
add_main_functions(BufrAscatObsBuilder)
