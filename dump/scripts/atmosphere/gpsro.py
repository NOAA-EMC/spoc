#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions, map_path
from gpsro_obs_builder import BaseGpsroBufrObsBuilder

MAPPING_PATH = map_path("./gpsro.yaml")

# ----------------------------------------------------------------------
# Concrete implementation
# ----------------------------------------------------------------------


class GpsroBufrObsBuilder(BaseGpsroBufrObsBuilder):
    """Supply the YAML files to the base class."""

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))


add_main_functions(GpsroBufrObsBuilder)
