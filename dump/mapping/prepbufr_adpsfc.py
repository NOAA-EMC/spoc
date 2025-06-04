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

    def _make_description():
        description = super()._make_description()

        description.add_variables([
            {
                'name': "MetaData/dateTime",
                'source': 'timestamp',
                'units': "seconds since 1970-01-01T00:00:00Z",
                'longName': "Observation Time"
            },
            {
                'name': 'MetaData/sequenceNumber',
                'source': 'sequenceNumber',
                'units': '1',
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
        container = bufr.Parser(input_path, mapping_path).parse(comm)

        self.log.debug(f'container list (original): {container.list()}')

        self.log.debug(f'Do DateTime calculation')
        otmct = container.get('obsTimeMinusCycleTime')
        otmct_paths = container.get_paths('obsTimeMinusCycleTime')
        otmct2 = np.array(otmct)

        self._add_timestamp(container, self._get_reference_time(input_path))

        self.log.debug(f'Make an array of 0s for MetaData/sequenceNumber')
        sequenceNum = np.zeros(lon.shape, dtype=np.int32)
        self.log.debug(f' sequenceNummin/max =  {sequenceNum.min()} {sequenceNum.max()}')

        self.log.debug(f'Add variables to container')
        container.add('sequenceNumber', sequenceNum, lon_paths)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')

        return container


add_main_functions(AdpsfcPrepbufrObsBuilder)
