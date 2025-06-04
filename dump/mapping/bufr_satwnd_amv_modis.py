#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions, map_path
from bufr_satwnd_amv_obs_builder import SatWndAmvObsBuilder


MAPPING_PATH = map_path('bufr_satwnd_amv_modis.yaml')
FIND_QI = 1


class SatWndAmvModisObsBuilder(SatWndAmvObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        container = super().make_obs(comm, input_path)

        # Add new/derived data into container
        for cat in container.all_sub_categories():
            self._add_quality_info_and_gen_app(FIND_QI, container, cat)

        return container

    def _make_description(self):
        description = super()._make_description()
        self._add_quality_info_and_gen_app_descriptions(description)

        return description

    def _get_obs_type(self, swcm, chan_freq):
        obstype = swcm.copy()

        # Use numpy vectorized operations
        obstype = np.where(swcm >= 4, 259, obstype)  # WVCA/DL
        obstype = np.where(swcm == 3, 258, obstype)  # WVCT
        obstype = np.where(swcm == 1, 257, obstype)  # IRLW

        if not np.any(np.isin(obstype, [257, 258, 259])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype


# Add main functions create_obs_file and create_obs_group
add_main_functions(SatWndAmvModisObsBuilder)
