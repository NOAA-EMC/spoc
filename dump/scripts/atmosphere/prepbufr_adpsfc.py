#!/usr/bin/env python3
import calendar
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.obs_builder import add_main_functions
from prepbufr_obs_builder import PrepbufrObsBuilder, map_path


MAPPING_PATH = map_path('prepbufr_adpsfc.yaml')


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
            }
        ])
        return description

    def make_obs(self, comm, input_path):
        """
        Create the ioda adpsfc prepbufr observations:
        - reads values
        - adds sequenceNum

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

        self.log.debug(f'Make an array of 0s for MetaData/sequenceNumber')
        sequenceNum = np.zeros(dhr.shape, dtype=np.int32)
        self.log.debug(f' sequenceNummin/max =  {sequenceNum.min()} {sequenceNum.max()}')

        self.log.debug(f'Do tsen and tv calculation')
        tpc = container.get('temperatureEventCode')
        tob = container.get('airTemperatureObsValue')
        tob_paths = container.get_paths('airTemperatureObsValue')
        tsen = np.full(tob.shape[0], tob.fill_value)
        tsen = np.where(((tpc >= 1) & (tpc < 8)), tob, tsen)
        tvo = np.full(tob.shape[0], tob.fill_value)
        tvo = np.where((tpc == 8), tob, tvo)

        self.log.debug(f'Do tsen and tv QM calculations')
        tobqm = container.get('airTemperatureQualityMarker')
        tsenqm = np.full(tobqm.shape[0], tobqm.fill_value)
        tsenqm = np.where(((tpc >= 1) & (tpc < 8)), tobqm, tsenqm)
        tvoqm = np.full(tobqm.shape[0], tobqm.fill_value)
        tvoqm = np.where((tpc == 8), tobqm, tvoqm)

        self.log.debug(f'Do tsen and tv ObsError calculations')
        toboe = container.get('airTemperatureObsError')
        tsenoe = np.full(toboe.shape[0], toboe.fill_value)
        tsenoe = np.where(((tpc >= 1) & (tpc < 8)), toboe, tsenoe)
        tvooe = np.full(toboe.shape[0], toboe.fill_value)
        tvooe = np.where((tpc == 8), toboe, tvooe)

        self.log.debug(f'Update variables in container')
        container.replace('airTemperatureObsValue', tsen)
        container.replace('airTemperatureQualityMarker', tsenqm)
        container.replace('airTemperatureObsError', tsenoe)
        container.replace('virtualTemperatureObsValue', tvo)
        container.replace('virtualTemperatureQualityMarker', tvoqm)
        container.replace('virtualTemperatureObsError', tvooe)

        self.log.debug(f'Add variables to container')
        container.add('sequenceNumber', sequenceNum, dhr_paths)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')

        return container


add_main_functions(AdpsfcPrepbufrObsBuilder)
