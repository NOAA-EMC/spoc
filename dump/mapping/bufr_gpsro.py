#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions, map_path
from bufr_gpsro_obs_builder import BaseGpsroBufrObsBuilder

LOCATION_PROFILE = map_path("./bufr_gpsro_latitude.yaml")
HEIGHT_PROFILE   = map_path("./bufr_gpsro_height.yaml")

# ----------------------------------------------------------------------
# Concrete implementation
# ----------------------------------------------------------------------
class GpsroBufrObsBuilder(BaseGpsroBufrObsBuilder):
    """Supply the two YAML files to the base class."""

    def __init__(self):
        map_dict = {
            "loc_profile":   str(LOCATION_PROFILE),
            "height_profile": str(HEIGHT_PROFILE),
        }
        super().__init__(map_dict, log_name=os.path.basename(__file__))


add_main_functions(GpsroBufrObsBuilder)
