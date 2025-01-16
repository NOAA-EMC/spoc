#!/usr/bin/env python3

import sys
from b2ibase.util import parse_arguments
from b2ibase.config import Config
from b2ibase.data_container import DataContainer
from b2ibase.b2i import B2I
from b2ibase.log import B2ILogger


class TesacData(DataContainer):
    def process_data(self):
        self.convert_temp_to_celsius("waterTemperature")

        temp = self.container.get("variables/waterTemperature")
        temp_min = -10.0
        temp_max = 50.0
        temp_mask = (temp > temp_min) & (temp < temp_max)

        saln = self.container.get("variables/salinity")
        saln_min = 0.0
        saln_max = 45.0
        saln_mask = (saln > saln_min) & (saln < saln_max)

        station_id = self.container.get("variables/stationID")
        id_mask = [True if id.isdigit() and id != 0 else False for id in station_id]
        mask = id_mask & temp_mask & saln_mask
        self.filter(mask)

        saln_error = 0.01
        temp_error = 0.02

        self.add_preqc_var("waterTemperature")
        self.add_preqc_var("salinity")
        self.add_error_var("waterTemperature", temp_error)
        self.add_error_var("salinity", saln_error)
        self.add_seq_num()


if __name__ == '__main__':

    script_name, config_file, log_file, test_file = parse_arguments()
    log_to_console = True
    logger = B2ILogger(script_name, log_to_console, log_file)

    config = Config(config_file, logger)
    data = TesacData(logger)
    b2i = B2I(config, data, logger)
    b2i.run()
    if test_file:
        result = b2i.test(test_file)
        sys.exit(result)
