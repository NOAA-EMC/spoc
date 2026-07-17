#!/usr/bin/env python3
import os
#import numpy as np
#import numpy.ma as ma

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path
#from retrieval_ozone_obs_builder import BufrOzoneObsBuilder

MAPPING_PATH = map_path('retrieval_lghtng_vaisala.yaml')


class BufrlghtngObsBuilder(ObsBuilder):
    """
    Class for building observations from ompst8 BUFR data.

    This class extends `ObsBuilder` to include specific logic for processing
    Level-2 retrived total ozone data from OMPS nadir mapper

    :param mapping_path: Path to the mapping file.
    :type mapping_path: str
    """

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))


# Add main functions create_obs_file or create_obs_group
add_main_functions(BufrlghtngObsBuilder)
