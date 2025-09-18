#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions, map_path
from bufr_satwnd_amv_obs_builder import SatWndAmvObsBuilder


MAPPING_PATH = map_path('bufr_satwnd_amv_leogeo.yaml')


class SatWndAmvLeogeoObsBuilder(SatWndAmvObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def _get_obs_type(self, swcm, chan_freq):
        obstype = swcm.copy()

        # Use numpy vectorized operations
        obstype = np.where(swcm == 1, 255, obstype)  # IRLW

        if not np.any(np.isin(obstype, [255])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype


add_main_functions(SatWndAmvLeogeoObsBuilder)
