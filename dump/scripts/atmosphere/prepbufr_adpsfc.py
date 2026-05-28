#!/usr/bin/env python3
import calendar
import os
import numpy as np
import numpy.ma as ma
import yaml

import bufr
from bufr.obs_builder import add_main_functions
from prepbufr_obs_builder import PrepbufrObsBuilder, map_path


MAPPING_PATH = map_path('prepbufr_adpsfc.yaml')

# Fixed number of temperature event levels handled by this builder.
NUM_T_EVENTS = 5

# - If ObsType/virtualTemperature is in the encoder variables, use Tv if available otherwise Tdry.
#   This option mimics the GSI's TSENSIBLE=False default
# - If ObsType/virtualTemperature not in encoder variables, always use Tdry
#   This is what we want to do long-term

# obs types 181, 187 (land stations) always use Tdry to match GSI behavior
TSENSIBLE_EXCEPTION_TYPES = [181, 187]


def _check_include_tv(yaml_path):
    """Check if virtualTemperature should be included based on encoder variables in YAML."""
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    encoder_vars = config.get('encoder', {}).get('variables', [])
    return any(v.get('name') == 'ObsType/virtualTemperature' for v in encoder_vars)

class AdpsfcPrepbufrObsBuilder(PrepbufrObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def _make_description(self):
        description = super()._make_description()

        variables = [
            {
                'name': 'MetaData/sequenceNumber',
                'source': 'sequenceNumber',
                'longName': 'Sequence Number (Obs Subtype)',
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

    def make_obs(self, comm, input_path):
        """
        Create the ioda adpsfc prepbufr observations:
        - reads values
        - adds sequenceNum
        - adds ObsSubType
        - extracts Tdry from event stack (peels back through events if top is Tv)

        Parameters
        ----------
        comm: object
                The communicator object (e.g., MPI)
        input_path: str
                The input bufr file
        """

        container = super().make_obs(comm, input_path)

        # Get container from mapping file first
        self.log.info('Get container from bufr')
        container = bufr.Parser(input_path, MAPPING_PATH).parse(comm)

        self.log.debug(f'container list (original): {container.list()}')

        self.log.debug(f'Do DateTime calculation')
        dhr = container.get('obsTimeMinusCycleTime')
        dhr_paths = container.get_paths('obsTimeMinusCycleTime')
        dhr2 = np.array(dhr)
        self._replace_timestamp(container, self._get_reference_time(input_path))

        self.log.debug(f'Make an array of 0s for MetaData/sequenceNumber and ObsSubType')
        sequenceNum = np.zeros(dhr.shape, dtype=np.int32)
        self.log.debug(f' sequenceNum min/max =  {sequenceNum.min()} {sequenceNum.max()}')

        include_tv = _check_include_tv(MAPPING_PATH)
        self.log.debug(f'Extract temperature from event stack (include_tv={include_tv})')

        # get record-specific data
        toboe = container.get('airTemperatureObsError')
        obs_type = container.get('observationType')

        # get event-specific data
        tpc_events = []
        tob_events = []
        tqm_events = []
        for i in range(1, NUM_T_EVENTS + 1):
            tpc_events.append(container.get(f'temperatureEventCode{i}'))
            tob_events.append(container.get(f'temperatureOb{i}'))
            tqm_events.append(container.get(f'temperatureQM{i}'))

        # get paths for adding new variables
        tob_paths = container.get_paths('temperatureOb1')

        # create arrays with fill_value (matching develop branch pattern)
        n_obs = tob_events[0].shape[0]
        tsen = np.full(n_obs, tob_events[0].fill_value)
        tsenqm = np.full(n_obs, tqm_events[0].fill_value)
        tsenoe = np.full(n_obs, toboe.fill_value)
        tvo = np.full(n_obs, tob_events[0].fill_value)
        tvoqm = np.full(n_obs, tqm_events[0].fill_value)
        tvooe = np.full(n_obs, toboe.fill_value)

        # loop through obs
        for idx in range(n_obs):
            use_tv = include_tv and (int(obs_type[idx]) not in TSENSIBLE_EXCEPTION_TYPES)

            # look back through events for desired T field
            for ev in range(NUM_T_EVENTS):
                tpc_val = tpc_events[ev][idx]
                tob_val = tob_events[ev][idx]
                tqm_val = tqm_events[ev][idx]

                if ma.is_masked(tpc_val) or ma.is_masked(tob_val):
                    continue

                # select desired obs type, if present 
                if tpc_val == 8 and use_tv:
                    # use Tv if available
                    tvo[idx] = tob_val
                    if not ma.is_masked(tqm_val):
                        tvoqm[idx] = tqm_val
                    if not ma.is_masked(toboe[idx]):
                        tvooe[idx] = toboe[idx]
                    break
                elif (tpc_val >= 1) and (tpc_val < 8):
                    # Save Tdry
                    tsen[idx] = tob_val
                    if not ma.is_masked(tqm_val):
                        tsenqm[idx] = tqm_val
                    if not ma.is_masked(toboe[idx]):
                        tsenoe[idx] = toboe[idx]
                    break

        self.log.debug(f'Update variables in container')
        container.add('airTemperatureObsValue', tsen, tob_paths)
        container.add('airTemperatureQualityMarker', tsenqm, tob_paths)
        container.replace('airTemperatureObsError', tsenoe)

        if include_tv:
            container.add('virtualTemperatureObsValue', tvo, tob_paths)
            container.add('virtualTemperatureQualityMarker', tvoqm, tob_paths)
            container.add('virtualTemperatureObsError', tvooe, tob_paths)

        self.log.debug(f'Add variables to container')
        container.add('sequenceNumber', sequenceNum, dhr_paths)
        container.add('obsSubType', sequenceNum, dhr_paths)

        self.log.debug(f'container list (updated): {container.list()}')

        return container


add_main_functions(AdpsfcPrepbufrObsBuilder)

