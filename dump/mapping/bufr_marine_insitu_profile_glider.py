#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions
from bufr_marine_insitu_obs_builder import MarineInsituObsBuilder, map_path

MAPPING_PATH = map_path('bufr_marine_insitu_profile_glider.yaml')


class MarineInsituProfileGliderObsBuilder(MarineInsituObsBuilder):
    def __init__(self, config=None):
        super().__init__(MAPPING_PATH, config=config, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
        container = super().make_obs(comm, input_path)

        temp = container.get("waterTemperature")
        temp_mask = (temp > -10.0) & (temp < 50.0)

        saln = container.get("salinity")
        saln_mask = (saln > 0.0) & (saln < 45.0)

        id = container.get("stationID")
        id_mask = (id >= 68900) & (id <= 68999) | \
                  (id >= 1800000) & (id <= 1809999) | \
                  (id >= 2800000) & (id <= 2809999) | \
                  (id >= 3800000) & (id <= 3809999) | \
                  (id >= 4800000) & (id <= 4809999) | \
                  (id >= 5800000) & (id <= 5809999) | \
                  (id >= 6800000) & (id <= 6809999) | \
                  (id >= 7800000) & (id <= 7809999)

        container.apply_mask(id_mask & temp_mask & saln_mask)

        self._add_preqc_var(container, "waterTemperature")
        self._add_preqc_var(container, "salinity")
        self._add_error_var(container, "waterTemperature", error=0.02)
        self._add_error_var(container, "salinity", error=0.01)
        self._add_seq_num(container)

        return container

add_main_functions(MarineInsituProfileGliderObsBuilder)
