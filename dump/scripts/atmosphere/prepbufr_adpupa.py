#!/usr/bin/env python3
import os
import numpy as np
import time
import calendar
from datetime import datetime

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path
from prepbufr_obs_builder import PrepbufrObsBuilder

MAPPING_PATH = map_path('prepbufr_adpupa.yaml')


class AdpupaPrepbufrObsBuilder(PrepbufrObsBuilder):
    """
    A builder class to generate ADPUPA observations from ADPUPA prepBUFR input.
    """

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def compute_conditional_array(self, source_array, condition_mask):
        """
        Compute an array where values from source_array are retained
        if condition_mask is True, else fill_value is used.
        """

        result = np.full(source_array.shape, source_array.fill_value)
        result[condition_mask] = source_array[condition_mask]
        return result

    def make_obs(self, comm, input_path):

        # Get container from mapping file first
        self.log.info(f'Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'container list (original): {container.list()}')

        self.log.debug(f'Perform DateTime calculation')
        hrdr = container.get('obsTimeMinusCycleTime')
        self._replace_timestamp(container, self._get_reference_time(input_path))

        self.log.debug(f'Perform stationPressure, stationPressureQM, and stationPressureError calculations')
        cat = container.get('prepbufrDataLevelCategory')
        pob = container.get('pressure')
        pqm = container.get('pressureQualityMarker')
        poe = container.get('pressureError')

        station_pressure = self.compute_conditional_array(pob, cat == 0)
        station_pressureQM = self.compute_conditional_array(pqm, cat == 0)
        station_pressureError = self.compute_conditional_array(poe, cat == 0)

        self.log.debug(f'Perform airTemperature, airTemperatureQM, and airTemperatureError calculations')
        tpc = container.get('temperatureEventCode')
        tob = container.get('airTemperature')
        tobqm = container.get('airTemperatureQualityMarker')
        toboe = container.get('airTemperatureError')

        air_temperature = self.compute_conditional_array(tob, (tpc >= 1) & (tpc < 8))
        air_temperatureQM = self.compute_conditional_array(tobqm, (tpc >= 1) & (tpc < 8))
        air_temperatureError = self.compute_conditional_array(toboe, (tpc >= 1) & (tpc < 8))

        self.log.debug(f'Perform virtualTemperature, virtualTemperatureQM, and virtualTemperatureError calculations')
        virtual_temperature = self.compute_conditional_array(tob, tpc == 8)
        virtual_temperatureQM = self.compute_conditional_array(tobqm, tpc == 8)
        virtual_temperatureError = self.compute_conditional_array(toboe, tpc == 8)

        self.log.debug(f'Update variables into container')
        container.replace('airTemperature', air_temperature)
        container.replace('virtualTemperature', virtual_temperature)
        container.replace('airTemperatureQualityMarker', air_temperatureQM)
        container.replace('airTemperatureError', air_temperatureError)
        container.replace('virtualTemperatureQualityMarker', virtual_temperatureQM)
        container.replace('virtualTemperatureError', virtual_temperatureError)

        self.log.debug(f'Add new/derived variables into container')
        ydr_paths = container.get_paths('latitude')
        container.add('stationPressure', station_pressure, ydr_paths)
        container.add('stationPressureQualityMarker', station_pressureQM, ydr_paths)
        container.add('stationPressureError', station_pressureError, ydr_paths)

        self.log.debug(f'container list (updated): {container.list()}')

        return container

    def _make_description(self):
        description = super()._make_description()

        description.add_variables([
            {
                'name': 'ObsValue/stationPressure',
                'source': 'stationPressure',
                'units': 'Pa',
                'longName': 'Station Pressure',
            },
            {
                'name': 'QualityMarker/stationPressure',
                'source': 'stationPressureQualityMarker',
                'units': '',
                'longName': 'Station Pressure Quality Marker',
            },
            {
                'name': 'ObsError/stationPressure',
                'source': 'stationPressureError',
                'units': 'Pa',
                'longName': 'Station Pressure Error',
            },
        ])

        return description


# Add main functions create_obs_file or create_obs_group
add_main_functions(AdpupaPrepbufrObsBuilder)
