#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions
from bufr_marine_insitu_obs_builder import MarineInsituObsBuilder, map_path

MAPPING_PATH = map_path('bufr_marine_insitu_surface_drifter.yaml')


class MarineInsituSurfaceDrifterObsBuilder(MarineInsituObsBuilder):
    def __init__(self, config=None):
        super().__init__(MAPPING_PATH, config=config, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
        container = super().make_obs(comm, input_path)

        temp = container.get("seaSurfaceTemperature")
        temp_mask = (temp > -10.0) & (temp < 50.0)

        buoy_type = container.get("buoyType")

        # Separate Drifter profiles from dbuoy tank
        # buoy_type:
        # 1 - Standard Lagrangian drifter (Global Drifter Programme)
        # 4 - Ice drifter
        # 5 - SVPG Standard Lagrangian drifter with GPS
        buoy_mask = np.isin(buoy_type,  [1, 4, 5])

        print("BBBBBBBBBBBBBBBBBBB")
        print(buoy_mask)
        print(temp_mask.size)
        print(buoy_type.size)
        print(buoy_mask.size)
        print(np.count_nonzero(buoy_mask))
        print(np.count_nonzero(temp_mask))
        print("BBBBBBBBBBBBBBBBBBB")
        container.apply_mask(buoy_mask & temp_mask)

        self._add_preqc_var(container, "seaSurfaceTemperature")
        self._add_preqc_var(container, "salinity")
        self._add_error_var(container, "seaSurfaceTemperature", error=0.24)
        self._add_error_var(container, "salinity", error=0.01)

        return container

add_main_functions(MarineInsituSurfaceDrifterObsBuilder)
