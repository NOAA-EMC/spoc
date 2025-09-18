#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import ObsBuilder, add_dummy_variable


class BufrOzoneObsBuilder(ObsBuilder):
    """
    Observations builder for BUFR Ozone data.

    Inherits from ObsBuilder and provides specific logic for ozone observation processing.
    """

    def __init__(self, mapping_path, log_name=os.path.basename(__file__)):
        super().__init__(mapping_path, log_name=log_name)

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

            self._add_pressures(container, cat)

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
                'name': 'MetaData/pressure',
                'source': 'pressure',
                'units': 'Pa',
                'longName': 'Pressure',
            }])

    def _add_pressures(self, container, cat):
        """
        Add pressure metadata to the container for a given category.

        :param container: Observation container.
        :type container: object
        :param cat: Category identifier (tuple or similar).
        :type cat: object
        """
        # Add new MetaData variables: MetaData/pressure
        satId = container.get('satelliteId', cat)
        if not satId.size:
            self.log.warning(f'category {cat[0]} does not exist in input file')
            add_dummy_variable(container, 'pressure', cat, 'latitude')
            return

        latitude = container.get('latitude', cat)
        paths = container.get_paths('latitude', cat)

        pressure = np.full_like(latitude, 0)
        container.add('pressure', pressure, paths, cat)
