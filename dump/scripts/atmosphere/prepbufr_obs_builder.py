#!/usr/bin/env python3

import os
import re
import numpy as np
import numpy.ma as ma
from pathlib import Path

from datetime import datetime

import bufr
from bufr.obs_builder import ObsBuilder


def map_path(map_file_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, map_file_name)


class PrepbufrObsBuilder(ObsBuilder):
    def __init__(self, mapping_path, log_name=os.path.basename(__file__)):
        super().__init__(mapping_path, log_name=log_name)

    def _get_reference_time(self, input_path) -> np.datetime64:
        path_components = Path(input_path).parts

        dump_regex = r'\w+\.(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})'
        test_regex = r'(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})(?P<hour>\d{2})'

        for idx, component in enumerate(reversed(path_components)):
            dump_match = re.match(dump_regex, component)
            test_match = re.match(test_regex, component)
            if dump_match:
                if idx == len(path_components) - 1:
                    continue

                ref_time = datetime(year=int(dump_match.group('year')),
                                    month=int(dump_match.group('month')),
                                    day=int(dump_match.group('day')),
                                    hour=int(path_components[-1*(idx+1) + 1]))
                break
            elif test_match:
                ref_time = datetime(year=int(test_match.group('year')),
                                    month=int(test_match.group('month')),
                                    day=int(test_match.group('day')),
                                    hour=int(test_match.group('hour')))
                break
        else:
            print(f'Reference date not found in path.')
            ref_time = datetime(year=2020, month=1, day=1)

        return np.datetime64(ref_time)

    def _compute_datetime(self, cycleTimeSinceEpoch, dhr):
        """
        Compute dateTime using the cycleTimeSinceEpoch and Observation Time
            minus Cycle Time

        Parameters:
            cycleTimeSinceEpoch: Time of cycle in Epoch Time
            dhr: Observation Time Minus Cycle Time

        Returns:
            Masked array of dateTime values
        """

        int64_fill_value = np.int64(0)

        dateTime = np.zeros(dhr.shape, dtype=np.int64)
        for i in range(len(dateTime)):
            if ma.is_masked(dhr[i]):
                continue
            else:
                dateTime[i] = np.int64(dhr[i]*3600) + cycleTimeSinceEpoch

        dateTime = ma.array(dateTime)
        dateTime = ma.masked_values(dateTime, int64_fill_value)

        return dateTime

    def _replace_timestamp(self, container: bufr.DataContainer, reference_time: np.datetime64) -> np.array:
        times = container.get('obsTimeMinusCycleTime')

        cycle_times = ma.masked_array(np.round(3600 * times).astype(np.int64), dtype='timedelta64[s]', mask=times.mask)
        timestamps = ma.masked_array(reference_time + cycle_times,
                                     mask=times.mask, dtype='datetime64[s]').astype('int64')

        container.replace('timestamp', timestamps.filled())
