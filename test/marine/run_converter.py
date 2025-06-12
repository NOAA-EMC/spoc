import os
import subprocess
from b2i_config import *
from b2i_tests import *
from marine_config import generate_test_case


def run_converter(test_case, converter_dir, input_dir, output_dir, config_dir):
    converter = os.path.join(converter_dir, test_case["converter"])
    input_path = os.path.join(input_dir, test_case["input"])
    output_path = os.path.join(output_dir, test_case["reference"])
    config_path = os.path.join(config_dir, test_case["config"])

    cmd = [converter, "--input", input_path, "--output", output_path, "--config", config_path]
    print(f"Running command: {' '.join(cmd)}", flush=True)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f'Success!')
        return True
    except subprocess.CalledProcessError as e:
        print(f'Error running  {cmd}')
        return False


if __name__ == "__main__":

    converter_dir = "/work/noaa/da/edwardg/spoc/dump/mapping/"
    input_dir = "/work/noaa/da/edwardg/spoc/test/marine_data/testdata"
    output_dir = "."
    config_dir = "/work/noaa/da/edwardg/spoc/test/marine_data/testconfig"

    # test_case = generate_test_case(CONFIG_TEST_DATA[0])
    # run_converter(test_case, converter_dir, input_dir, output_dir, config_dir)

    for test in CONFIG_TEST_DATA:
        test_case = generate_test_case(test)
        run_converter(test_case, converter_dir, input_dir, output_dir, config_dir)
