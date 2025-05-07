#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import ObsBuilder


def map_path(map_file_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, map_file_name)


class SatWndAmvObsBuilder(ObsBuilder):
    def __init__(self, mapping_path, log_name=os.path.basename(__file__)):
        super().__init__(mapping_path, log_name=log_name)

    def make_obs(self, comm, input_path):
        # Get container from mapping file first
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'container list (original): {container.list()}')
        self.log.debug(f'all_sub_categories =  {container.all_sub_categories()}')
        self.log.debug(f'category map =  {container.get_category_map()}')

        # Add new/derived data into container
        for cat in container.all_sub_categories():
            self.log.debug(f'category = {cat}')

            satId = container.get('satelliteId', cat)
            if not np.any(satId):
                self.log.warning(f'category {cat[0]} does not exist in input file')

            self._add_wind_obs(container, cat)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')
        self.log.debug('all_sub_categories {container.all_sub_categories()}')

        return container

    # Provide defualt implementations for methods from the ObsBuilder class
    def _make_description(self):
        description = super()._make_description()
        self._add_wind_descriptions(description)

        return description

    # Methods that are used to extend the export description
    def _add_wind_descriptions(self, description):
        description.add_variables([
            {
                'name': 'ObsType/windEastward',
                'source': 'obstype_uwind',
                'units': '1',
                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
            },
            {
                'name': 'ObsType/windNorthward',
                'source': 'obstype_vwind',
                'units': '1',
                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
            },
            {
                'name': 'ObsValue/windEastward',
                'source': 'windEastward',
                'units': 'm s-1',
                'longName': 'Eastward Wind Component',
            },
            {
                'name': 'ObsValue/windNorthward',
                'source': 'windNorthward',
                'units': 'm s-1',
                'longName': 'Northward Wind Component',
            }])

    def _add_quality_info_and_gen_app_descriptions(self, description):
        description.add_variables([
            {
                'name': 'MetaData/windGeneratingApplication',
                'source': 'windGeneratingApplication',
                'units': '1',
                'longName': 'Wind Generating Application',
            },
            {
                'name': 'MetaData/qiWithoutForecast',
                'source': 'qualityInformationWithoutForecast',
                'units': '1',
                'longName': 'Quality Information Without Forecast',
            }])

    # Methods that are used to extend the obs data container
    def _add_wind_obs(self, container, cat):
        # Add new variables: ObsType/windEastward & ObsType/windNorthward
        swcm = container.get('windComputationMethod', cat)
        chanfreq = container.get('sensorCentralFrequency', cat)

        if swcm.size == 0:
            self.log.warning(f'category {cat[0]} does not exist in input file')
            paths = container.get_paths('variables/windComputationMethod', cat)
            obstype = container.get('variables/windComputationMethod', cat)
            container.add('variables/obstype_uwind', obstype, paths, cat)
            container.add('variables/obstype_vwind', obstype, paths, cat)

            paths = container.get_paths('variables/windSpeed', cat)
            wob = container.get('variables/windSpeed', cat)
            container.add('variables/windEastward', wob, paths, cat)
            container.add('variables/windNorthward', wob, paths, cat)
            return

        # self.log.debug(f'swcm min/max = {swcm.min()} {swcm.max()}')
        self.log.debug('chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')

        obstype = self._get_obs_type(swcm, chanfreq)

        self.log.debug(f'obstype = {obstype}')
        self.log.debug(f'obstype min/max =  {obstype.min()} {obstype.max()}')

        paths = container.get_paths('windComputationMethod', cat)
        container.add('obstype_uwind', obstype, paths, cat)
        container.add('obstype_vwind', obstype, paths, cat)

        # Add new variables: ObsValue/windEastward & ObsValue/windNorthward
        wdir = container.get('windDirection', cat)
        wspd = container.get('windSpeed', cat)

        self.log.debug(f'wdir min/max = {wdir.min()} {wdir.max()}')
        self.log.debug(f'wspd min/max = {wspd.min()} {wspd.max()}')

        uob, vob = self._compute_wind_components(wdir, wspd)

        self.log.debug(f'uob min/max = {uob.min()} {uob.max()}')
        self.log.debug(f'vob min/max = {vob.min()} {vob.max()}')

        paths = container.get_paths('windSpeed', cat)
        container.add('windEastward', uob, paths, cat)
        container.add('windNorthward', vob, paths, cat)

    def _add_quality_info_and_gen_app(self, findQi, container, cat):
        # Add new variables: MetaData/windGeneratingApplication and qiWithoutForecast
        gnap2D = container.get('generatingApplication', cat)
        pccf2D = container.get('qualityInformation', cat)
        satId = container.get('satelliteId', cat)

        if not satId.size:
            paths = container.get_paths('windComputationMethod', cat)
            dummy = container.get('windSpeed', cat)
            container.add('windGeneratingApplication', dummy, paths, cat)
            container.add('qualityInformationWithoutForecast', dummy, paths, cat)
            return

        gnap, qifn = self._get_quality_info_and_gen_app(findQi, gnap2D, pccf2D)

        self.log.debug(f'gnap min/max = {gnap.min()} {gnap.max()}')
        self.log.debug(f'qifn min/max = {qifn.min()} {qifn.max()}')

        paths = container.get_paths('windComputationMethod', cat)
        container.add('windGeneratingApplication', gnap, paths, cat)
        container.add('qualityInformationWithoutForecast', qifn, paths, cat)

    def _get_obs_type(self, swcm, chan_freq=0):
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

        raise NotImplementedError('Method _get_obs_type must be implemented in derived classes')

    # Private methods
    def _compute_wind_components(self, wdir, wspd):
        """
        Compute the U and V wind components from wind direction and wind speed.

        Parameters:
            wdir (array-like): Wind direction in degrees (meteorological convention: 0° = North, 90° = East).
            wspd (array-like): Wind speed.

        Returns:
            tuple: U and V wind components as numpy arrays with dtype float32.
        """
        wdir_rad = np.radians(wdir)  # Convert degrees to radians
        u = -wspd * np.sin(wdir_rad)
        v = -wspd * np.cos(wdir_rad)

        return u.astype(np.float32), v.astype(np.float32)

    def _get_quality_info_and_gen_app(self, findQi, gnap2D, pccf2D):
        # For NOAA VIIRS data, qi w/o forecast (qifn) is packaged in same
        # vector of qi with ga = 5 (EUMETSAT QI without forecast). Must
        # conduct a search and extract the correct vector for gnap and qi

        # 1. Find dimension-sizes of ga and qi (should be the same!)
        gDim1, gDim2 = np.shape(gnap2D)
        qDim1, qDim2 = np.shape(pccf2D)
        self.log.info('Generating Application and Quality Information SEARCH')
        self.log.debug(f'Dimension size of GNAP ({gDim1},{gDim2})')
        self.log.debug(f'Dimension size of PCCF ({qDim1},{qDim2})')

        # 2. Initialize gnap and qifn as None, and search for dimension of
        #    ga with values of 5. If the same column exists for qi, assign
        #    gnap to ga[:,i] and qifn to qi[:,i], else raise warning that no
        #    appropriate GNAP/PCCF combination was found
        gnap = None
        qifn = None
        for i in range(gDim2):
            if np.unique(gnap2D[:, i].squeeze()) == findQi:
                if i <= qDim2:
                    self.log.info(f'GNAP/PCCF found for column {i}')
                    gnap = gnap2D[:, i].squeeze()
                    qifn = pccf2D[:, i].squeeze()
                else:
                    self.log.info(f'ERROR: GNAP column {i} outside of PCCF dimension {qDim2}')
        if (gnap is None) & (qifn is None):
            raise ValueError(f'GNAP == {findQI} NOT FOUND OR OUT OF PCCF DIMENSION-RANGE, WILL FAIL!')
        # If EE is needed, key search on np.unique(gnap2D[:,i].squeeze()) == 7 instead
        # NOTE: Make sure to return np.float32 or np.int32 types as appropriate!!!
        return gnap.astype(np.int32), qifn.astype(np.int32)
