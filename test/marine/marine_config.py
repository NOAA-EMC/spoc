import os
import yaml
from b2i_config import *
from b2i_tests import *


# utility to create yaml configuration files
# for bufr to marine ioda converters

def create_test_config_files(b2i_config, test_dir):
    testconfig_dir = os.path.join(test_dir, "testconfig")
    testdata_dir = os.path.join(test_dir, "testdata")
    testresult_dir = os.path.join(test_dir, "testresults")

    for test in CONFIG_TEST_DATA:
        i = test.data_type
        testconfig_path = os.path.join(testconfig_dir, b2i_config_filenames[i])
        print(f'----->>>>  create_test_config_files: {testconfig_path}')
        # print(f'----->>>>  creating {test}')
        # print(f'----->>>>  creating {test.data_format}')

        b2i_config.create_config_file(
            test.data_format,
            test.subsets,
            test.data_type,
            test.data_description,
            cycle_type,
            cycle_datetime,
            testdata_dir,
            testresult_dir,
            OCEAN_BASIN_FILE,
            testconfig_path
        )

def generate_test_case(test):
    i = test.data_type
    if test.data_type in marine_profile_instruments:
        descriptor = "profile_" + test.data_type
    elif test.data_type in marine_surface_instruments:
        descriptor = "surface_" + test.data_type
    else:
        print(f"Error: unknown data_type {test.data_type}")

    return {
        "name": b2i_test_names[i],
        "converter": b2i_converters[i],
        "input": bufr_filename(cycle_datetime, cycle_type, cycle, test.data_format),
        "reference": ioda_filename(cycle_type, cycle, descriptor, cycle_datetime),
        "config": b2i_config_filenames[i]
    }


def generate_marine_test_config_file(config_filename, test_dir):
    test_cases = []

    for test in CONFIG_TEST_DATA:
        test_cases.append(generate_test_case(test))

    converter_dir = "/work/noaa/da/edwardg/spoc/dump/mapping/"

    # Define the test suite structure
    test_suite = {
        "name": "marine",
        "converter_dir": converter_dir,
        "test_data_dir": test_dir,
        "tests": test_cases
    }

    # Define the YAML structure
    yaml_data = {
        "test_suites": [test_suite]
    }

    # Write to YAML file
    with open(config_filename, "w") as yaml_file:
        yaml.dump(yaml_data, yaml_file, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":

    b2i_config = Bufr2iodaConfig()

    test_dir = "/work/noaa/da/edwardg/spoc/test/marine_data"
    create_test_config_files(b2i_config, test_dir)

    config_filename = 'marine_tests.yaml'
    generate_marine_test_config_file(config_filename, test_dir)
