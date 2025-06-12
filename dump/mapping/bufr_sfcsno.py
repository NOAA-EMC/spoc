#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions


MAPPING_PATH = map_path('bufr_sfcsno.yaml')


class BufrSfcsnoObsBuilder(ObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'container list (original): {container.list()}')

        # Add new/derived data into container
        sogr = container.get('groundState')
        snod = container.get('totalSnowDepth')
        snod[(sogr <= 11.0) & snod.mask] = 0.0
        snod[(sogr == 15.0) & snod.mask] = 0.0
        snod.mask = (snod < 0.0) | snod.mask
        container.replace('totalSnowDepth', snod)
        snod_upd = container.get('totalSnowDepth')

        container.apply_mask(~snod.mask)

        return container


add_main_functions(BufrSfcsnoObsBuilder)
