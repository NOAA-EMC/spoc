#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions
from ..config_path import config_path


MAPPING_PATH = config_path("atmosphere", "radiance_amsr2.yaml")


class BufrAmsr2ObsBuilder(ObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))


add_main_functions(BufrAmsr2ObsBuilder)
