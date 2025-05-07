#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions
from bufr_satwnd_amv_obs_builder import SatWndAmvObsBuilder, map_path


MAPPING_PATH = map_path('bufr_satwnd_amv_ahi.yaml')
FIND_QI = 102


class SatWndAmvAhiObsBuilder(SatWndAmvObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
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
        obstype = np.where(swcm == 5, 250, obstype)  # WVCA/DL
        obstype = np.where(swcm == 3, 250, obstype)  # WVCT
        obstype = np.where(swcm == 2, 242, obstype)  # VIS
        obstype = np.where(swcm == 1, 252, obstype)  # IRLW

        if not np.any(np.isin(obstype, [242, 250, 252])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype


add_main_functions(SatWndAmvAhiObsBuilder)
