#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions
from bufr_marine_insitu_obs_builder import MarineInsituObsBuilder, map_path

MAPPING_PATH = map_path('bufr_marine_insitu_surface_dbuoyb_drifter.yaml')


'''
buoy types for drifters:
------------------------
00      Unspecified drifting buoy
01      Standard Lagrangian drifter (Global Drifter Programme)
02      Standard FGGE type drifting buoy (non-Lagrangian meteorological drifting buoy)
03      Wind measuring FGGE type drifting buoy (non-Lagrangian meteorological drifting buoy)
04      Ice drifter
05      SVPG Standard Lagrangian drifter with GPS (BUFR)
06      SVP-HR drifter with high-resolution temperature or thermistor string (BUFR)
10      ALACE (Autonomous Lagrangian Circulation Explorer)
11      MARVOR (MARine VORtical profiler)
12      RAFOS (Ranging and Fixing of Sound)
13      PROVOR (Profiling float with Argos)
14      SOLO (Swimbladder-Operated Lagrangian Oscillating)
15      APEX (Autonomous Profiling Explorer)
'''

drifter_buoy_types = [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15] 


class MarineInsituSurfaceDrifterObsBuilder(MarineInsituObsBuilder):
    def __init__(self, config=None):
        super().__init__(MAPPING_PATH, config=config, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
        container = super().make_obs(comm, input_path)

        temp = container.get("seaSurfaceTemperature")
        temp_mask = (temp > -10.0) & (temp < 50.0)

        buoy_type = container.get("buoyType")
        rpid = container.get("stationID")

        # rpid = stationID: string array (e.g., 'A8xxx')
        # buoy_type: int array (e.g., 1, 2, 3), etc.
        drifter_mask = np.isin(buoy_type, drifter_buoy_types, assume_unique=True)

        # Optional: Add RPID check for drifter patterns (e.g., starts with 'A8')
        rpid_drifter_mask = np.array([isinstance(r, str) and r.startswith('A8') for r in rpid.filled('')])
        drifter_mask = drifter_mask | rpid_drifter_mask

        # Handle masked (missing) BUYT values
        # If BUYT is masked, assume not a drifter unless RPID suggests otherwise
        drifter_mask = np.where(buoy_type.mask, rpid_drifter_mask, drifter_mask)

        # print("BBBBBBBBBBBBBBBBBBB")
        # print(buoy_mask)
        # print(temp_mask.size)
        # print(buoy_type.size)
        # print(buoy_mask.size)
        # print(np.count_nonzero(buoy_mask))
        # print(np.count_nonzero(temp_mask))
        # print("BBBBBBBBBBBBBBBBBBB")
        # container.apply_mask(buoy_mask & temp_mask)

        container.apply_mask(drifter_mask & temp_mask)

        self._add_preqc_var(container, "seaSurfaceTemperature")
        self._add_error_var(container, "seaSurfaceTemperature", error=0.24)

        return container

add_main_functions(MarineInsituSurfaceDrifterObsBuilder)
