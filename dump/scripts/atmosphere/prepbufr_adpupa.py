#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma
import time
import calendar
import yaml
from datetime import datetime

import bufr
from bufr.obs_builder import ObsBuilder, add_main_functions, map_path
from prepbufr_obs_builder import PrepbufrObsBuilder

MAPPING_PATH = map_path('prepbufr_adpupa.yaml')
NUM_T_EVENTS = 5


def _check_include_tv(yaml_path):
    """Check if virtualTemperature should be included based on encoder variables in YAML."""
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    encoder_vars = config.get('encoder', {}).get('variables', [])
    return any(v.get('name') == 'ObsType/virtualTemperature' for v in encoder_vars)


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

        include_tv = _check_include_tv(MAPPING_PATH)
        self.log.debug(f'Extract temperature from event stack (include_tv={include_tv})')

        toboe = container.get('airTemperatureError')
        tpc_events = []
        tob_events = []
        tqm_events = []
        for i in range(1, NUM_T_EVENTS + 1):
            tpc_events.append(container.get(f'temperatureEventCode{i}'))
            tob_events.append(container.get(f'temperatureOb{i}'))
            tqm_events.append(container.get(f'temperatureQM{i}'))

        n_obs = tob_events[0].shape[0]
        air_temperature = np.full(n_obs, tob_events[0].fill_value)
        air_temperatureQM = np.full(n_obs, tqm_events[0].fill_value)
        air_temperatureError = np.full(n_obs, toboe.fill_value)
        derived_temperature_event_code = np.full(n_obs, tpc_events[0].fill_value)

        if include_tv:
            virtual_temperature = np.full(n_obs, tob_events[0].fill_value)
            virtual_temperatureQM = np.full(n_obs, tqm_events[0].fill_value)
            virtual_temperatureError = np.full(n_obs, toboe.fill_value)

        for idx in range(n_obs):
            selected_tdry = None
            selected_tv = None

            for ev in range(NUM_T_EVENTS):
                tpc_val = tpc_events[ev][idx]
                tob_val = tob_events[ev][idx]
                tqm_val = tqm_events[ev][idx]

                if ma.is_masked(tpc_val) or ma.is_masked(tob_val):
                    continue

                if selected_tdry is None and (tpc_val >= 1) and (tpc_val < 8):
                    selected_tdry = (tpc_val, tob_val, tqm_val)
                    if not include_tv:
                        break
                if include_tv and selected_tv is None and (tpc_val == 8):
                    selected_tv = (tpc_val, tob_val, tqm_val)
                    if selected_tdry is not None:
                        break

            if selected_tdry is not None:
                tpc_val, tob_val, tqm_val = selected_tdry
                air_temperature[idx] = tob_val
                if not ma.is_masked(tqm_val):
                    air_temperatureQM[idx] = tqm_val
                if not ma.is_masked(toboe[idx]):
                    air_temperatureError[idx] = toboe[idx]

            if include_tv:
                selected_output = selected_tv if selected_tv is not None else selected_tdry
                if selected_output is not None:
                    tpc_val, tob_val, tqm_val = selected_output
                    virtual_temperature[idx] = tob_val
                    if not ma.is_masked(tqm_val):
                        virtual_temperatureQM[idx] = tqm_val
                    if not ma.is_masked(toboe[idx]):
                        virtual_temperatureError[idx] = toboe[idx]
                    # With one metadata field, temperatureEventCode tracks the selected
                    # virtual-temperature output when enabled (Tv with Tdry fallback).
                    derived_temperature_event_code[idx] = tpc_val
            elif selected_tdry is not None:
                derived_temperature_event_code[idx] = selected_tdry[0]

        self.log.debug(f'Update variables into container')
        container.replace('airTemperature', air_temperature)
        container.replace('airTemperatureQualityMarker', air_temperatureQM)
        container.replace('airTemperatureError', air_temperatureError)
        container.replace('temperatureEventCode', derived_temperature_event_code)
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

        if _check_include_tv(MAPPING_PATH):
            variables.append({
                'name': 'ObsSubType/virtualTemperature',
                'source': 'obsSubType',
                'longName': 'Observation SubType',
            })

        description.add_variables(variables)

        return description


# Add main functions create_obs_file or create_obs_group
add_main_functions(AdpupaPrepbufrObsBuilder)
