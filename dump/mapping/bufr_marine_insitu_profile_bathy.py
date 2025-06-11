#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions
from bufr_marine_insitu_obs_builder import MarineInsituObsBuilder, map_path

MAPPING_PATH = map_path('bufr_marine_insitu_profile_bathy.yaml')


class MarineInsituProfileBathyObsBuilder(MarineInsituObsBuilder):
    def __init__(self, config=None):
        super().__init__(MAPPING_PATH, config=config, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
        container = super().make_obs(comm, input_path)

        temp = container.get("waterTemperature")
        container.apply_mask((temp > -10.0) & (temp < 50.0))

        self._add_preqc_var(container, "waterTemperature")
        self._add_error_var(container, "waterTemperature", error=0.24)
        self._add_seq_num(container)

        return container

add_main_functions(MarineInsituProfileBathyObsBuilder)
