#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path


MAPPING_PATH = map_path('radiance_gmi.yaml')


class BufrAtmsObsBuilder(ObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))


add_main_functions(BufrAtmsObsBuilder)
