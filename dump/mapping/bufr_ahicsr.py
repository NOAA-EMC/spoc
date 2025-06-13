#!/usr/bin/env python3
import os
import numpy as np

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path

MAPPING_PATH = map_path('bufr_ahicsr.yaml')


class BufrAhicsrObsBuilder(ObsBuilder):

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def compute_sensor_channel_number(self, nlocs, nchannels):
        """
        Creates a 2D array of sensor channel numbers with shape (nlocs, nchannels),
        where each row is [1, 2, ..., nchannels].

        Args:
            nlocs (int): Number of observation locations.
            nchannels (int): Number of channels (default is 12).

        Returns:
            np.ndarray: 2D array with shape (nlocs, nchannels).
        """
        return np.tile(np.arange(7, nchannels + 1, dtype=np.int32), (nlocs, 1))

    def make_obs(self, comm, input_path):
        # Get container from mapping file
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'Container list (original): {container.list()}')
        self.log.debug(f'All_sub_categories =  {container.all_sub_categories()}')
        self.log.debug(f'Category map = {container.get_category_map()}')

        # Add a new derived variable, sensorChannelNumber into container
        for cat in container.all_sub_categories():
            self.log.debug(f'category: {cat}')

            satId = container.get('satelliteId', cat)
            sccf_paths = container.get_paths('sensorCentralFrequency', cat)
            nlocs = satId.shape[0]  # Number of locations
            nchannels = 16  # Number of channels 
            # Add Ten Channels from 7 to 16 
            sensor_channel_number = self.compute_sensor_channel_number(nlocs, nchannels)
            self.log.debug(f'Adding derived variable: sensorChannelNumber for {nlocs} locations')
            container.add('sensorChannelNumber', sensor_channel_number, sccf_paths, cat)

            if not np.any(satId):
                self.log.warning(f'Category {cat[0]} does not exist in input file')
                continue  # Skip invalid category

        # Final container state
        self.log.debug(f'Container list (updated): {container.list()}')

        return container

    def _make_description(self):
        description = super()._make_description()

        description.add_variables([
            {
                'name': 'MetaData/sensorChannelNumber',
                'source': 'sensorChannelNumber',
                'units': '',
                'longName': 'Sensor Channel Number',
            },
        ])

        return description


add_main_functions(BufrAhicsrObsBuilder)
