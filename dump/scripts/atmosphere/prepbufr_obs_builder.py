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

        if num_events == 0:
            return tsen, tsenqm, tsenoe, tvo, tvoqm, tvooe

        # Stack the event-stack levels into (num_events, n_obs) arrays so the
        # search over levels is a handful of NumPy ops instead of a Python
        # loop over every observation.
        tpc_stack = ma.stack(tpc_events[:num_events])
        tob_stack = ma.stack(tob_events[:num_events])
        tqm_stack = ma.stack(tqm_events[:num_events])

        valid = ~ma.getmaskarray(tpc_stack) & ~ma.getmaskarray(tob_stack)
        tpc_filled = ma.filled(tpc_stack, -1)

        is_tv_code = valid & (tpc_filled == 8) & use_tv_arr[np.newaxis, :]
        is_sen_code = valid & (tpc_filled >= 1) & (tpc_filled < 8)
        match = is_tv_code | is_sen_code

        # First matching event level per observation (top of stack wins);
        # argmax returns the first True, and 0 (harmless) when none match.
        has_match = match.any(axis=0)
        first_ev = np.argmax(match, axis=0)
        obs_idx = np.arange(n_obs)

        sel_is_tv = has_match & is_tv_code[first_ev, obs_idx]
        sel_is_sen = has_match & ~sel_is_tv

        tob_sel = ma.filled(tob_stack, tob_events[0].fill_value)[first_ev, obs_idx]
        tqm_sel = ma.filled(tqm_stack, tqm_events[0].fill_value)[first_ev, obs_idx]
        tqm_valid = ~ma.getmaskarray(tqm_stack)[first_ev, obs_idx]
        toboe_valid = ~ma.getmaskarray(toboe)
        toboe_filled = ma.filled(toboe, toboe.fill_value)

        tsen[sel_is_sen] = tob_sel[sel_is_sen]
        tsenqm[sel_is_sen & tqm_valid] = tqm_sel[sel_is_sen & tqm_valid]
        tsenoe[sel_is_sen & toboe_valid] = toboe_filled[sel_is_sen & toboe_valid]

        tvo[sel_is_tv] = tob_sel[sel_is_tv]
        tvoqm[sel_is_tv & tqm_valid] = tqm_sel[sel_is_tv & tqm_valid]
        tvooe[sel_is_tv & toboe_valid] = toboe_filled[sel_is_tv & toboe_valid]

        return tsen, tsenqm, tsenoe, tvo, tvoqm, tvooe
