#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions, map_path
from bufr_gpsro_obs_builder import BaseGpsroBufrObsBuilder

MAPPING_PATH = map_path("./bufr_gpsro.yaml")

# ----------------------------------------------------------------------
# Concrete implementation
# ----------------------------------------------------------------------


class GpsroBufrObsBuilder(BaseGpsroBufrObsBuilder):
    print("NICKE 1")
    """Supply the two YAML files to the base class."""

    def __init__(self):
        print("NICKE x2")
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))


add_main_functions(GpsroBufrObsBuilder)
