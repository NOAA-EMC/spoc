#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import add_main_functions, map_path, add_dummy_variable
from bufr_satwnd_amv_obs_builder import SatWndAmvObsBuilder


MAPPING_PATH = map_path('bufr_satwnd_amv_avhrr.yaml')


class SatWndAmvAvhrrObsBuilder(SatWndAmvObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
        container = super().make_obs(comm, input_path)

        # Add new/derived data into container
        for cat in container.all_sub_categories():
            self._add_avhrr_quality_info_and_gen_app(container, cat)

        return container

    def _make_description(self):
        description = super()._make_description()
        self._add_quality_info_and_gen_app_descriptions(description)

        return description

    def _get_obs_type(self, swcm, chan_freq):
        obstype = swcm.copy()

        # Use numpy vectorized operations
        obstype = np.where(swcm == 1, 244, obstype)  # IRLW

        if not np.any(np.isin(obstype, [244])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype.astype(np.int32)

    def _add_avhrr_quality_info_and_gen_app(self, container, cat):
        # Add new variables: MetaData/windGeneratingApplication and qualityInformationWithoutForecast
        gnap2D = container.get('generatingApplication', cat)
        pccf2D = container.get('qualityInformation', cat)
        satId = container.get('satelliteId', cat)

        if not satId.size:

            dummy_mappings = [
                ('windGeneratingApplication', 'windComputationMethod'),
                ('qualityInformationWithoutForecast', 'windSpeed')
            ]
            for target_var, source_var in dummy_mappings:
                add_dummy_variable(container, target_var, cat, source_var)

            return

        gnap, qifn = self._get_avhrr_quality_info_and_gen_app(gnap2D, pccf2D, satId)

        self.log.debug(f'gnap min/max = {gnap.min()} {gnap.max()}')
        self.log.debug(f'qifn min/max = {qifn.min()} {qifn.max()}')

        paths = container.get_paths('windComputationMethod', cat)
        container.add('windGeneratingApplication', gnap, paths, cat)
        container.add('qualityInformationWithoutForecast', qifn, paths, cat)

    def _get_avhrr_quality_info_and_gen_app(self, gnap2D, pccf2D, satID):
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
            # If findQI is not found anywhere in gnap2D, set findQI and findEE to 1 and 4, respectively
            if not np.any(np.isin(gnap2D, [findQI])):
                self.log.debug(
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
            if np.unique(gnap2D[:, i]) == findQI:
                if i < qDim2:
                    self.log.info(f'GNAP/PCCF found for column {i}')
                    gnap = gnap2D[:, i].copy()
                    qifn = pccf2D[:, i].copy()
                else:
                    self.log.info(f'ERROR: GNAP column {i} outside of PCCF dimension {qDim2}')
        if (gnap is None) and (qifn is None):
            raise ValueError(f'GNAP == {findQI} NOT FOUND OR OUT OF PCCF DIMENSION-RANGE, WILL FAIL!')
        # If EE is needed, key search on np.unique(gnap2D[:,i].squeeze()) == findEE instead
        return gnap, qifn


# Add main functions create_obs_file and create_obs_group
add_main_functions(SatWndAmvAvhrrObsBuilder)
