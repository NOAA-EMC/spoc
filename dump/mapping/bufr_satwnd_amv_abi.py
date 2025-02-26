#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions

from bufr_satwnd_amv_obs_builder import SatWndAmvObsBuilder


class SatWndAmvAbiObsBuilder(SatWndAmvObsBuilder):
    def __init__(self, input_path, mapping_path):
        super().__init__(input_path, mapping_path, log_name=os.path.basename(__file__))

    # Use the SatWndAmvObsBuilder implementations for the methods make_description and make_obs

    def _get_obs_type(self, swcm, chanfreq):
        """
        Determine the observation type based on `swcm` and `chanfreq`.

        Parameters:
            swcm (array-like): Satellite derived wind calculation method.
            chanfreq (array-like): Satellite channel center frequency (Hz).

        Returns:
            numpy.ndarray: Observation type array.

        Raises:
            ValueError: If any `obstype` is unassigned.
        """

        obstype = swcm.copy()

        # Use numpy vectorized operations
        obstype = np.where(swcm == 5, 247, obstype)  # WVCA/DL
        obstype = np.where(swcm == 3, 246, obstype)  # WVCT
        obstype = np.where(swcm == 2, 251, obstype)  # VIS
        obstype = np.where(swcm == 1, 245, obstype)  # IRLW

        condition = np.logical_and(swcm == 1, chanfreq >= 5e13)  # IRSW
        obstype = np.where(condition, 240, obstype)

        if not np.any(np.isin(obstype, [247, 246, 251, 245, 240])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype


# Add main functions create_obs_file and create_obs_group
add_main_functions(SatWndAmvAbiObsBuilder)
