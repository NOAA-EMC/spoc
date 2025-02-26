#!/usr/bin/env python3

import os
import numpy as np

import bufr
from bufr.obs_builder import ObsBuilder


class SatWndAmvObsBuilder(ObsBuilder):
    def __init__(self, input_path, mapping_path, log_name=os.path.basename(__file__)):
        super().__init__(input_path, mapping_path, log_name=log_name)

    # Provide defualt implementations for methods from the ObsBuilder class
    def _make_description(self):
        description = bufr.encoders.Description(self.mapping_path)
        self.add_wind_descriptions(description)
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

            self.add_wind_obs(container, cat)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')
        self.log.debug('all_sub_categories {container.all_sub_categories()}')

        return container


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
                'name': 'qualityInformationWithoutForecast',
                'source': 'qualityInformationWithoutForecast',
                'units': '1',
                'longName': 'Quality Information Without Forecast',
            }])


    # Methods that are used to extend the obs data container
    def _add_wind_obs(self, container, cat):
        # Add new variables: ObsType/windEastward & ObsType/windNorthward
        swcm = container.get('windComputationMethod', cat)
        chanfreq = container.get('sensorCentralFrequency', cat)

        self.log.debug(f'swcm min/max = {swcm.min()} {swcm.max()}')
        self.log.debug('chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')

        obstype = self._get_obs_type(swcm)

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

    def _add_quality_info_and_gen_app(self, container, cat):
        # Add new variables: MetaData/windGeneratingApplication and qualityInformationWithoutForecast
        gnap2D = container.get('generatingApplication', cat)
        pccf2D = container.get('qualityInformation', cat)
        satId = container.get('satelliteId', cat)

        gnap, qifn = self._get_quality_info_and_gen_app(gnap2D, pccf2D, satId)

        self.log.debug(f'gnap min/max = {gnap.min()} {gnap.max()}')
        self.log.debug(f'qifn min/max = {qifn.min()} {qifn.max()}')

        paths = container.get_paths('windComputationMethod', cat)
        container.add('windGeneratingApplication', gnap, paths, cat)
        container.add('qualityInformationWithoutForecast', qifn, paths, cat)


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


    def _get_obs_type(self, swcm):
        raise NotImplementedError("Method '_get_obs_type' must be implemented in a subclass")

    def _get_quality_info_and_gen_app(self, gnap2D, pccf2D, satId):
        raise NotImplementedError("Method '_get_quality_info_and_gen_app' must be implemented in a subclass")
