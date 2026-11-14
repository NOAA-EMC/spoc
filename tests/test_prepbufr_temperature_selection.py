import importlib
import sys
import types
from pathlib import Path

import numpy as np
import numpy.ma as ma
import yaml


def _import_prepbufr_obs_builder(monkeypatch):
    fake_bufr = types.ModuleType('bufr')
    fake_obs_builder = types.ModuleType('bufr.obs_builder')

    class FakeObsBuilder:
        def __init__(self, *args, **kwargs):
            pass

    class FakeDataContainer:
        pass

    fake_obs_builder.ObsBuilder = FakeObsBuilder
    fake_bufr.obs_builder = fake_obs_builder
    fake_bufr.DataContainer = FakeDataContainer

    monkeypatch.setitem(sys.modules, 'bufr', fake_bufr)
    monkeypatch.setitem(sys.modules, 'bufr.obs_builder', fake_obs_builder)

    scripts_path = Path(__file__).resolve().parents[1] / 'dump' / 'scripts' / 'atmosphere'
    monkeypatch.syspath_prepend(str(scripts_path))

    if 'prepbufr_obs_builder' in sys.modules:
        del sys.modules['prepbufr_obs_builder']

    return importlib.import_module('prepbufr_obs_builder')


def _build_temperature_inputs():
    tpc_events = [
        ma.array([8, 8, 2], fill_value=99),
        ma.array([1, 1, 1], fill_value=99),
    ]
    tob_events = [
        ma.array([301.0, 311.0, 290.0], fill_value=-9999.0),
        ma.array([300.0, 310.0, 289.0], fill_value=-9999.0),
    ]
    tqm_events = [
        ma.array([1, 2, 3], fill_value=255),
        ma.array([4, 5, 6], fill_value=255),
    ]
    toboe = ma.array([0.5, 0.6, 0.7], fill_value=-9999.0)
    use_tv = np.array([True, False, True])
    return tpc_events, tob_events, tqm_events, toboe, use_tv


def test_check_legacy_tv_selection_default_and_enabled(tmp_path, monkeypatch):
    prepbufr_obs_builder = _import_prepbufr_obs_builder(monkeypatch)

    cfg_default = {'bufr': {'subsets': ['ADPSFC']}}
    cfg_enabled = {'bufr': {'subsets': ['ADPSFC'], 'legacy_tv_over_tdry': True}}

    default_path = tmp_path / 'default.yaml'
    enabled_path = tmp_path / 'enabled.yaml'
    default_path.write_text(yaml.safe_dump(cfg_default), encoding='utf-8')
    enabled_path.write_text(yaml.safe_dump(cfg_enabled), encoding='utf-8')

    assert prepbufr_obs_builder.check_legacy_tv_selection(default_path) is False
    assert prepbufr_obs_builder.check_legacy_tv_selection(enabled_path) is True


def test_select_temperature_events_default_behavior(monkeypatch):
    prepbufr_obs_builder = _import_prepbufr_obs_builder(monkeypatch)
    tpc_events, tob_events, tqm_events, toboe, use_tv = _build_temperature_inputs()

    tsen, tsenqm, tsenoe, tvo, tvoqm, tvooe = prepbufr_obs_builder.PrepbufrObsBuilder._select_temperature_events(
        None, tpc_events, tob_events, tqm_events, toboe, use_tv, 2, legacy_tv_over_tdry=False)

    assert np.allclose(tsen, [300.0, 310.0, 290.0])
    assert np.array_equal(tsenqm, [4, 5, 3])
    assert np.allclose(tsenoe, [0.5, 0.6, 0.7])
    assert np.isclose(tvo[0], 301.0)
    assert np.isclose(tvoqm[0], 1)
    assert np.isclose(tvooe[0], 0.5)
    assert np.isclose(tvo[1], tob_events[0].fill_value)
    assert np.isclose(tvo[2], tob_events[0].fill_value)


def test_select_temperature_events_legacy_tv_over_tdry(monkeypatch):
    prepbufr_obs_builder = _import_prepbufr_obs_builder(monkeypatch)
    tpc_events, tob_events, tqm_events, toboe, use_tv = _build_temperature_inputs()

    tsen, tsenqm, tsenoe, tvo, tvoqm, tvooe = prepbufr_obs_builder.PrepbufrObsBuilder._select_temperature_events(
        None, tpc_events, tob_events, tqm_events, toboe, use_tv, 2, legacy_tv_over_tdry=True)

    assert np.isclose(tsen[0], tob_events[0].fill_value)
    assert np.isclose(tsenqm[0], tqm_events[0].fill_value)
    assert np.isclose(tsenoe[0], toboe.fill_value)
    assert np.allclose(tsen[1:], [310.0, 290.0])
    assert np.array_equal(tsenqm[1:], [5, 3])
    assert np.allclose(tsenoe[1:], [0.6, 0.7])
    assert np.isclose(tvo[0], 301.0)
    assert np.isclose(tvoqm[0], 1)
    assert np.isclose(tvooe[0], 0.5)
