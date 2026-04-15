#!/usr/bin/env python3
import calendar
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import add_main_functions
from prepbufr_obs_builder import PrepbufrObsBuilder, map_path


MAPPING_PATH = map_path('prepbufr_adpsfc.yaml')

# Number of temperature event levels read from YAML
NUM_T_EVENTS = 5

# Flag to mimi GSI's TSENSIBLE option. 
# True: Use sensible/dry temperature (Tdry) by searching through event stack
# False: Use Tv if available (TPC=8), otherwise Tdry (TPC 1-7), 
#        for testing against GSI only.
TSENSIBLE = True


class AdpsfcPrepbufrObsBuilder(PrepbufrObsBuilder):
    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def _make_description(self):
        description = super()._make_description()

        description.add_variables([
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
                'name': 'ObsSubType/virtualTemperature',
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
        ])

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

        self.log.debug(f'Extract temperature from event stack (tsensible={TSENSIBLE})')

        # Get temperature data from all event levels
        tpc_events = []
        tob_events = []
        tqm_events = []
        for i in range(1, NUM_T_EVENTS + 1):
            tpc_events.append(container.get(f'temperatureEventCode{i}'))
            tob_events.append(container.get(f'temperatureOb{i}'))
            tqm_events.append(container.get(f'temperatureQM{i}'))

        # Get ObsError (from T__BACKG, not event-specific)
        toboe = container.get('airTemperatureObsError')

        # Get fill values from each variable type
        tob_fill = tob_events[0].fill_value
        tqm_fill = tqm_events[0].fill_value
        toe_fill = toboe.fill_value

        # Initialize output arrays with appropriate fill values
        n_obs = tob_events[0].shape[0]
        tsen = np.full(n_obs, tob_fill)
        tsenqm = np.full(n_obs, tqm_fill)
        tsenoe = np.full(n_obs, toe_fill)
        tvo = np.full(n_obs, tob_fill)
        tvoqm = np.full(n_obs, tqm_fill)
        tvooe = np.full(n_obs, toe_fill)

        # Search through events for each observation
        for idx in range(n_obs):
            for ev in range(NUM_T_EVENTS):
                tpc_val = tpc_events[ev][idx]
                tob_val = tob_events[ev][idx]
                tqm_val = tqm_events[ev][idx]

                # Skip if obs masked/missing
                if ma.is_masked(tpc_val) or ma.is_masked(tob_val):
                    continue

                if tpc_val == 8 and ( not TSENSIBLE ):
                    # use Tv if available
                    tvo[idx] = tob_val
                    tvoqm[idx] = tqm_val
                    tvooe[idx] = toboe[idx]
                    break
                elif (tpc_val >= 1) and (tpc_val < 8):
                    # Save Tdry
                    tsen[idx] = tob_val
                    tsenqm[idx] = tqm_val
                    tsenoe[idx] = toboe[idx]
                    break

        self.log.debug(f'Update variables in container')
        container.replace('airTemperatureObsValue', tsen)
        container.replace('airTemperatureQualityMarker', tsenqm)
        container.replace('airTemperatureObsError', tsenoe)
        container.replace('virtualTemperatureObsValue', tvo)
        container.replace('virtualTemperatureQualityMarker', tvoqm)
        container.replace('virtualTemperatureObsError', tvooe)

        self.log.debug(f'Add variables to container')
        container.add('sequenceNumber', sequenceNum, dhr_paths)
        container.add('obsSubType', sequenceNum, dhr_paths)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')

        return container


add_main_functions(AdpsfcPrepbufrObsBuilder)

