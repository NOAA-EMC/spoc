#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma
import time
import calendar
from datetime import datetime

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path
from prepbufr_obs_builder import PrepbufrObsBuilder, check_include_tv

MAPPING_PATH = map_path('prepbufr_adpupa.yaml')
NUM_T_EVENTS = 5


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

        self.log.debug(f'Make an array of 0s for ObsSubType')
        obsSubType = np.zeros(hrdr.shape, dtype=np.int32)
        self.log.debug(f' obsSubType min/max =  {obsSubType.min()} {obsSubType.max()}')

        self.log.debug(f'Perform stationPressure, stationPressureQM, and stationPressureError calculations')
        cat = container.get('prepbufrDataLevelCategory')
        pob = container.get('pressure')
        pqm = container.get('pressureQualityMarker')
        poe = container.get('pressureError')

        station_pressure = self.compute_conditional_array(pob, cat == 0)
        station_pressureQM = self.compute_conditional_array(pqm, cat == 0)
        station_pressureError = self.compute_conditional_array(poe, cat == 0)

        include_tv = check_include_tv(MAPPING_PATH)
        self.log.debug(f'Extract temperature from event stack (include_tv={include_tv})')

        toboe = container.get('airTemperatureError')
        tpc_events = []
        tob_events = []
        tqm_events = []
        for i in range(1, NUM_T_EVENTS + 1):
            tpc_events.append(container.get(f'temperatureEventCode{i}'))
            tob_events.append(container.get(f'temperatureOb{i}'))
            tqm_events.append(container.get(f'temperatureQM{i}'))

        air_temperature, air_temperatureQM, air_temperatureError, \
            virtual_temperature, virtual_temperatureQM, virtual_temperatureError = \
            self._select_temperature_events(
                tpc_events, tob_events, tqm_events, toboe, include_tv, NUM_T_EVENTS)

        self.log.debug(f'Update variables into container')
        container.replace('airTemperature', air_temperature)
        container.replace('airTemperatureQualityMarker', air_temperatureQM)
        container.replace('airTemperatureError', air_temperatureError)
        if include_tv:
            container.replace('virtualTemperature', virtual_temperature)
            container.replace('virtualTemperatureQualityMarker', virtual_temperatureQM)
            container.replace('virtualTemperatureError', virtual_temperatureError)

        self.log.debug(f'Add new/derived variables into container')
        ydr_paths = container.get_paths('latitude')
        container.add('stationPressure', station_pressure, ydr_paths)
        container.add('stationPressureQualityMarker', station_pressureQM, ydr_paths)
        container.add('stationPressureError', station_pressureError, ydr_paths)
        container.add('obsSubType', obsSubType, ydr_paths)

        self.log.debug(f'container list (updated): {container.list()}')

        return container

    def _make_description(self):
        description = super()._make_description()

        variables = [
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
            {
                'name': 'ObsSubType/stationPressure',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            },
            {
                'name': 'ObsSubType/airTemperature',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            },
            {
                'name': 'ObsSubType/specificHumidity',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            },
            {
                'name': 'ObsSubType/windEastward',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            },
            {
                'name': 'ObsSubType/windNorthward',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            }
        ]

        if check_include_tv(MAPPING_PATH):
            variables.append({
                'name': 'ObsSubType/virtualTemperature',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            })

        description.add_variables(variables)

        return description


# Add main functions create_obs_file or create_obs_group
add_main_functions(AdpupaPrepbufrObsBuilder)
