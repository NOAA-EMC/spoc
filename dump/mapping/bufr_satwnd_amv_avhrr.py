#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions

from builders.satwnd_amv_obs_builder import SatWndAmvObsBuilder


class SatWndAmvAvhrrObsBuilder(SatWndAmvObsBuilder):
    def __init__(self, input_path, mapping_path):
        super().__init__(input_path, mapping_path, log_name=os.path.basename(__file__))

    # Override implementations for methods from the ObsBuilder class
    def _make_description(self):
        description = bufr.encoders.Description(self.mapping_path)

        self._add_wind_descriptions(description)
        self._add_quality_info_and_gen_app_descriptions(description)

        return description

    def _make_obs(self, comm):
        # Get container from mapping file first
        self.log.info('Get container from bufr')
        container = bufr.Parser(self.input_path, self.mapping_path).parse(comm)

        self.log.debug(f'container list (original): {container.list()}')
        self.log.debug(f'all_sub_categories =  {container.all_sub_categories()}')
        self.log.debug(f'category map =  {container.get_category_map()}')

        # Add new/derived data into container
        for cat in container.all_sub_categories():
            self.log.debug(f'category = {cat}')

            satId = container.get('satelliteId', cat)
            if not satId:
                self.log.warning(f'category {cat[0]} does not exist in input file')

            self._add_wind_obs(container, cat)
            self._add_gen_info_and_quality_info(container, cat)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')
        self.log.debug('all_sub_categories {container.all_sub_categories()}')

        return container

    # Override methods from SatWndAmvObsBuilder
    def _get_quality_info_and_gen_app(self, gnap2D, pccf2D, satID):
        # For METOP-A/B/C AVHRR data (satID 3,4,5), qi w/o forecast (qifn) is
        # packaged in same vector of qi with ga = 5 (QI without forecast), and EE
        # is packaged in same vector of qi with ga=7 (Estimated Error (EE) in m/s
        # converted to a percent confidence) shape (4,nobs).
        #
        # For NOAA-15/18/19 AVHRR data (satID 206,209,223), qi w/o forecast
        # (qifn) is packaged in same vector of qi with ga = 1 (EUMETSAT QI
        # without forecast), and EE is packaged in same vector of qi with ga=4
        # (Estimated Error (EE) in m/s converted to a percent confidence) shape
        # (4,nobs).
        #
        # Must conduct a search and extract the correct vector for gnap and qi
        # 0. Define the appropriate QI and EE search values, based on satID
        if np.all(np.isin(satID, [206, 209, 223])):  # NESDIS AVHRR set
            findQI = 1
            findEE = 4
        elif np.all(np.isin(satID, [3, 4, 5])):  # EUMETSAT AVHRR set
            findQI = 5
            findEE = 7
            # There is a catch: prior to 2023 AVHRR winds from EUMETSAT were formatted in the same
            # way as NESDIS AVHRR winds and both were passed through the NC005080 tank as a single
            # dataset. In that case, we need to actually set findQI=1 and findEE=4 here.
            # Let's do a preliminary check to see if any gnap2D values match findQI. If not, let's
            # automatically switch to findQI=1, findEE=4 and presume pre-2023 EUMETSAT AVHRR format
            if np.any(np.isin(gnap2D, [findQI])) == False:
                logging(comm, 'DEBUG',
                        f'NO GNAP VALUE OF {findQI} EXISTS FOR EUMETSAT AVHRR DATASET, PRESUMING PRE-2023 FORMATTING')
                findQI = 1
                findEE = 4
        else:
            self.log.debug(f'satID set not found (all satID values follow):')
            for sid in np.unique(satID):
                self.log.debug(f'satID: {sid}')
        self.log.debug(f'BTH: findQI={findQI}')
        # 1. Find dimension-sizes of ga and qi (should be the same!)
        gDim1, gDim2 = np.shape(gnap2D)
        qDim1, qDim2 = np.shape(pccf2D)
        self.log.info(f'Generating Application and Quality Information SEARCH:')
        self.log.debug(f'Dimension size of GNAP ({gDim1},{gDim2})')
        self.log.debug(f'Dimension size of PCCF ({qDim1},{qDim2})')
        # 2. Initialize gnap and qifn as None, and search for dimension of
        #    ga with values of findQI. If the same column exists for qi, assign
        #    gnap to ga[:,i] and qifn to qi[:,i], else raise warning that no
        #    appropriate GNAP/PCCF combination was found
        gnap = None
        qifn = None
        for i in range(gDim2):
            if np.unique(gnap2D[:, i].squeeze()) == findQI:
                if i <= qDim2:
                    self.log.info(f'GNAP/PCCF found for column {i}')
                    gnap = gnap2D[:, i].squeeze()
                    qifn = pccf2D[:, i].squeeze()
                else:
                    self.log.info(f'ERROR: GNAP column {i} outside of PCCF dimension {qDim2}')
        if (gnap is None) & (qifn is None):
            raise ValueError(f'GNAP == {findQI} NOT FOUND OR OUT OF PCCF DIMENSION-RANGE, WILL FAIL!')
        # If EE is needed, key search on np.unique(gnap2D[:,i].squeeze()) == findEE instead
        # NOTE: Make sure to return np.float32 or np.int32 types as appropriate!!!
        return gnap.astype(np.int32), qifn.astype(np.int32)

    def _get_obs_type(self, swcm):
        """
        Determine the observation type based on `swcm` and `chanfreq`.

        Parameters:
            swcm (array-like): Switch mode values.
            chanfreq (array-like): Channel frequency values (Hz).

        Returns:
            numpy.ndarray: Observation type array.

        Raises:
            ValueError: If any `obstype` is unassigned.
        """

        obstype = swcm.copy()

        # Use numpy vectorized operations
        obstype = np.where(swcm == 1, 244, obstype)  # IRLW

        if not np.any(np.isin(obstype, [244])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype.astype(np.int32)

# Add main functions create_obs_file and create_obs_group
add_main_functions(SatWndAmvAvhrrObsBuilder)
