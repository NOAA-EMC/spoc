#!/usr/bin/env python3
import os
from bufr.obs_builder import add_main_functions, map_path

from radiance_atovs_obs_builder import AtovsObsBuilder, NMFD, RARS


NMFD_MAPPING = map_path('radiance_amsua_1bamua.yaml')
RARS_MAPPING = map_path('radiance_amsua_esamua.yaml')


class AmsuaObsBuilder(AtovsObsBuilder):
    """
    ObsBuilder subclass for AMSU-A satellite data.

    Handles mapping, parsing, correction, and merging of AMSU-A 1B and ESA data
    using their respective mapping files.
    """

    def __init__(self):
        """
        Initialize the AmsuaObsBuilder.

        Sets up mapping dictionaries for 1B and ESA data types.
        """

        map_dict = {NMFD: NMFD_MAPPING,
                    RARS: RARS_MAPPING}

        super().__init__(map_dict, log_name=os.path.basename(__file__))


add_main_functions(AmsuaObsBuilder)
