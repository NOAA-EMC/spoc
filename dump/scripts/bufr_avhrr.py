#!/usr/bin/env python3
import os
import numpy as np

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path

MAPPING_PATH = map_path('bufr_avhrr.yaml')


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

            satId = container.get('satelliteId', cat)

            if not np.any(satId):
                self.log.warning(f'Category {cat[0]} does not exist in input file')
                continue  # Skip invalid category

        # Final container state
        self.log.debug(f'Container list (updated): {container.list()}')

        return container

    def _make_description(self):
        description = super()._make_description()
        return description


add_main_functions(BufrGsrasrObsBuilder)
