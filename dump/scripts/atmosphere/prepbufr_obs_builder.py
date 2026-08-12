#!/usr/bin/env python3

import os
import re
import numpy as np
import numpy.ma as ma
import yaml
from pathlib import Path

from datetime import datetime

import bufr
from bufr.obs_builder import ObsBuilder


def map_path(map_file_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, map_file_name)


def check_include_tv(yaml_path):
    """Check if virtualTemperature should be included based on encoder variables in YAML."""
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f) or {}
    encoder_vars = config.get('encoder', {}).get('variables', [])
    return any(v.get('name', '').startswith('ObsType/virtualTemperature') for v in encoder_vars)


class PrepbufrObsBuilder(ObsBuilder):
    def __init__(self, mapping_path, log_name=os.path.basename(__file__)):
        super().__init__(mapping_path, log_name=log_name)

    def _get_reference_time(self, input_path) -> np.datetime64:
        path_components = Path(input_path).parts
       """Extract date and hour from the directory path.  
          Looking for YYYYMMDDHH or YYYYMMDD/HH, with optional 
          "cycle." prefix."""

        ref_regex = re.compile(
            r'(?P<prefix>\w+\.)?'
            r'(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})'
            r'(?P<hour>\d{2})?$'
        )

        ref_time = None

        for idx, component in enumerate(reversed(path_components)):
            match = ref_regex.match(component)
            if not match:
                continue

            if match.group('hour') is not None:
                hour = int(match.group('hour'))
            else:
                # Date-only dir: get hour from next component.
                if idx == 0:
                    continue
                hour = int(path_components[-1 * (idx + 1) + 1])

            ref_time = datetime(
                year=int(match.group('year')),
                month=int(match.group('month')),
                day=int(match.group('day')),
                hour=hour
            )
            break

        if ref_time is None:
            self.log.error(f"Reference date not found in path: {input_path}")
            sys.exit(1)

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

        cycle_times = ma.masked_array(np.round(3600 * times).astype(np.int64),
                                      dtype='timedelta64[s]',
                                      mask=times.mask)

        timestamps = ma.masked_array(reference_time + cycle_times,
                                     mask=times.mask,
                                     fill_value=bufr.get_missing_value(np.int64),
                                     dtype='datetime64[s]').astype('int64')

        container.replace('timestamp', timestamps)

    def _select_temperature_events(self, tpc_events, tob_events, tqm_events, toboe, use_tv, num_events):
        """
        Select a single reported air temperature per observation from a
        stack of PREPBUFR temperature events, mirroring the GSI's Tsensible
        option: prefer the virtual-temperature event (temperatureEventCode
        == 8) if present and desired, otherwise fall back to the first
        sensible (Tdry) event (1 <= temperatureEventCode < 8). Exactly one
        of the sensible/virtual outputs is populated per observation -
        never both, since a Tv event is derived from an underlying Tdry
        event and both are otherwise present in the same stack.

        Parameters
        ----------
        tpc_events, tob_events, tqm_events: list of masked arrays
            Event-stack values (temperatureEventCode, temperatureOb,
            temperatureQM), one array per event level, ordered from the top
            of the stack.
        toboe: masked array
            Per-observation temperature obs error (not stacked by event).
        use_tv: bool or (n_obs,) bool array
            Whether virtual temperature should be preferred. Pass a
            per-observation array to exclude specific obs types (e.g. land
            stations) even when virtual temperature is enabled overall.
        num_events: int
            Number of event-stack levels to search.

        Returns
        -------
        tsen, tsenqm, tsenoe, tvo, tvoqm, tvooe: np.ndarray
            Sensible/virtual temperature, QM, and error arrays, each
            fill_value where not selected.
        """

        n_obs = tob_events[0].shape[0]
        use_tv_arr = np.broadcast_to(np.asarray(use_tv), (n_obs,))

        tsen = np.full(n_obs, tob_events[0].fill_value)
        tsenqm = np.full(n_obs, tqm_events[0].fill_value)
        tsenoe = np.full(n_obs, toboe.fill_value)
        tvo = np.full(n_obs, tob_events[0].fill_value)
        tvoqm = np.full(n_obs, tqm_events[0].fill_value)
        tvooe = np.full(n_obs, toboe.fill_value)

        for idx in range(n_obs):
            use_tv_idx = bool(use_tv_arr[idx])

            for ev in range(num_events):
                tpc_val = tpc_events[ev][idx]
                tob_val = tob_events[ev][idx]
                tqm_val = tqm_events[ev][idx]

                if ma.is_masked(tpc_val) or ma.is_masked(tob_val):
                    continue

                # select desired obs type, if present
                if tpc_val == 8 and use_tv_idx:
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

        return tsen, tsenqm, tsenoe, tvo, tvoqm, tvooe
