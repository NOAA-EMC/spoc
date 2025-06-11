import os
import json
from dataclasses import dataclass


OCEAN_BASIN_FILE = "/work/noaa/global/glopara/fix/gdas/soca/20240802/common/RECCAP2_region_masks_all_v20221025.nc"

cycle_type = "gdas"
cycle_datetime = '2019010700'
cycle = "00"


marine_profile_instruments = [
    "argo",
    "bathy",
    "glider",
    "tesac",
    "tropical",
    "xbtctd"
]
marine_surface_instruments = [
    "altkob",
    "cstgd",
    "drifter",
    "dbuoyb_drifter",
    "lcman",
    "shipsu",
    "trkob"
]

all_instruments = marine_profile_instruments + marine_surface_instruments


b2i_test_names = {}
for instrument in all_instruments:
    b2i_test_names[instrument] = "b2i_test_" + instrument


b2i_converters = {}
for instrument in marine_profile_instruments:
    b2i_converters[instrument] = "bufr_marine_insitu_profile_" + instrument + ".py"
for instrument in marine_surface_instruments:
    b2i_converters[instrument] = "bufr_marine_insitu_surface_" + instrument + ".py"

# print(json.dumps(b2i_converters, indent=4))



b2i_config_filenames = {}
for instrument in marine_profile_instruments:
    b2i_config_filenames[instrument] = "bufr2ioda_insitu_profile_" + instrument + "_" + cycle_datetime + ".yaml"
for instrument in marine_surface_instruments:
    b2i_config_filenames[instrument] = "bufr2ioda_insitu_surface_" + instrument + "_" + cycle_datetime + ".yaml"

# print(json.dumps(b2i_config_filenames, indent=4))


@dataclass
class ConfigTestData:
    data_format: str
    subsets: str
    data_type: str
    data_description: str

    # converter_filename: str
    # bufr_filename: str
    # ioda_filename: str
    # config_filename: str
    # surface_or_profile: str


# data_format, subsets, data_type
CONFIG_TEST_DATA = [
    ConfigTestData("subpfl", "SUBPFL", "argo", 
        '6-hrly in-situ ARGO profiles from subpfl: temperature and salinity'),
    ConfigTestData("bathy", "BATHY", "bathy", 
        '6-hrly in-situ profiles from BATHYthermal temperature'),
    ConfigTestData("subpfl", "SUBPFL", "glider", 
        '6-hrly in-situ Glider profiles from subpfl: temperature and salinity'),
    ConfigTestData("tesac", "TESAC", "tesac", 
        '6-hrly in-situ profiles from TESAC: temperature and salinity'),
    ConfigTestData("dbuoy", "DBUOY", "tropical", 
        '6-hrly in-situ tropical mooring profiles from dbuoy: temperature and salinity'),
    ConfigTestData("xbtctd", "XBTCTD", "xbtctd",  
        '6-hrly in-situ profiles from XBT/CTD: temperature and salinity'),
    ConfigTestData("altkob", "ALTKOB", "altkob", 
        '6-hrly in-situ surface obs from altkob: temperature and salinity'),
    ConfigTestData("cstgd", "CSTGD", "cstgd", 
        "6-hrly in-situ sea surface temperature obs from cstgd"),
    # ConfigTestData("dbuoy", "DBUOY", "drifter", 
        # "6-hrly in-situ Lagrangian drifter drogue profiles from dbuy: temperature"),
    # ConfigTestData("dbuoyb", "DBUOYB", "drifter", 
        # "6-hrly in-situ Lagrangian drifter drogue profiles from dbuoyb: temperature"),
    ConfigTestData("dbuoyb", "DBUOYB", "dbuoyb_drifter", 
        "6-hrly in-situ Lagrangian drifter drogue profiles from dbuoyb: temperature"),
    ConfigTestData("lcman", "LCMAN", "lcman", 
        '6-hrly in-situ surface temperature obs from LCMAN'),
    ConfigTestData("shipsu", "SHIPSU", "shipsu", 
        "6-hrly in-situ temperature obs from shipsu"),
    ConfigTestData("trkob", "TRACKOB", "trkob",  
        '6-hrly in-situ surface obs from TRACKOB: temperature and salinity')
]
