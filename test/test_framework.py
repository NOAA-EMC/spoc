import pytest
import yaml
import os
import tarfile
import requests
from pathlib import Path
import shutil
import subprocess
from run_compare import run_nccmp


def load_test_suites(yaml_file):
    # Check if the file exists
    if not os.path.exists(yaml_file):
        pytest.fail(f"YAML file '{yaml_file}' not found. Please provide a valid file using --test-config-file.")
    
    # Load the YAML file
    try:
        with open(yaml_file, 'r') as file:
            data = yaml.safe_load(file)
    except yaml.YAMLError as e:
        pytest.fail(f"Error parsing YAML file '{yaml_file}': {e}")
    except Exception as e:
        pytest.fail(f"Failed to read YAML file '{yaml_file}': {e}")

    # Validate the YAML structure
    if not isinstance(data, dict) or "test_suites" not in data:
        pytest.fail(f"Invalid YAML structure in '{yaml_file}'. Expected 'test_suites' key.")
    
    test_suites = data["test_suites"]
    if not test_suites:
        pytest.fail(f"No test suites found in '{yaml_file}'.")

    return test_suites


def download_and_extract_tarball(suite, downloads_dir):
    url = suite["url"]
    tarball = suite["tarball"]
    tarball_path = Path("test_suites") / tarball
    print(f"Downloading {tarball} from {url}")
    response = requests.get(f"{url}/{tarball}", stream=True)
    if response.status_code != 200:
        raise Exception(f"Failed to download {tarball}: HTTP {response.status_code}")
    tarball_path.parent.mkdir(exist_ok=True)
    with open(tarball_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    print(f"Extracting {tarball} to {downloads_dir}")
    with tarfile.open(tarball_path, "r:gz") as tar:
        tar.extractall(downloads_dir)
    tarball_path.unlink()  # Remove tarball after extraction


def create_subdir_symlink(symlink_path, subdir_path):
    # Check if symlink already exists
    if os.path.exists(symlink_path):
        if os.path.islink(symlink_path) and os.readlink(symlink_path) == subdir_path:
            print(f"Symlink already exists at {symlink_path} pointing to {subdir_path}.")
            return
        else:
            print(f"Error: {symlink_path} already exists but is not a symlink to {subdir_path}.")
            return

    if not os.path.exists(subdir_path):
        print(f"Error: {subdir_path} does not exist.")
    
    # Create the symlink
    try:
        os.symlink(subdir_path, symlink_path)
        print(f"Successfully created symlink at {symlink_path} pointing to {subdir_path}.")
    except OSError as e:
        print(f"Error creating symlink: {e}")


# Fixture to handle setup and teardown for each test suite
@pytest.fixture(scope="module")
def test_suite_setup(request):
    test_suite = request.param
    test_suite_dir = Path("test_suites") / test_suite["name"]

    # required subdir structure:
    # symlinks:
    testdata_dir = test_suite_dir / "testdata"
    testoutput_dir = test_suite_dir / "testoutput"
    testconfig_dir = test_suite_dir / "testconfig"  # optional
    # directories:
    testresults_dir = test_suite_dir / "testresults"
    downloads_dir = test_suite_dir / "downloads"    # if tar ball

    # if test_suite_dir exists, and has the right structure, use it.
    # otherwise, create it and place in it symlinks to user's data,
    # which is either in a given dir or in some dir unpacked from
    # a tar ball

    if test_suite_dir.exists():
        # check that it has the required subdirectories
        if not os.path.exists(testdata_dir):
            pytest.fail(f"Setup failed: {testdata_dir} does not exist.")
        if not os.path.exists(testoutput_dir):
            pytest.fail(f"Setup failed: {testoutput_dir} does not exist.")
        print(f"Using existing data in {test_suite_dir}")
    else:
        os.makedirs(test_suite_dir, exist_ok=True)

        if "test_data_dir" in test_suite:
            user_test_dir = Path(test_suite["test_data_dir"])
        else:
            os.makedirs(downloads_dir, exist_ok=True)
            download_and_extract_tarball(test_suite, downloads_dir)
            user_test_dir = downloads_dir

        if user_test_dir.exists():
            print(f"Using data in {user_test_dir}")
        else:
            pytest.fail(f"Setup failed: {user_test_dir} does not exist.")

        user_input_dir = user_test_dir / "testdata"
        user_reference_dir = user_test_dir / "testoutput"
        user_config_dir = user_test_dir / "testconfig"
        user_results_dir = user_test_dir / "testresults"

        if not os.path.exists(user_input_dir):
            pytest.fail(f"Setup failed: {user_input_dir} does not exist.")
        create_subdir_symlink(testdata_dir, user_input_dir)

        if not os.path.exists(user_reference_dir):
            pytest.fail(f"Setup failed: {user_reference_dir} does not exist.")
        create_subdir_symlink(testoutput_dir, user_reference_dir)

        if os.path.exists(user_config_dir):
            create_subdir_symlink(testconfig_dir, user_config_dir)

        if os.path.exists(user_results_dir):
            create_subdir_symlink(testresults_dir, user_results_dir)
    os.makedirs(testresults_dir, exist_ok=True)

    yield {
        "suite_name": test_suite["name"],
        "testdata_dir": testdata_dir,
        "testoutput_dir": testoutput_dir,
        "testconfig_dir": testconfig_dir,
        "testresults_dir": testresults_dir,
        "converter_dir": test_suite.get("converter_dir"),
        "tests": test_suite["tests"],
    }

    # Optional cleanup (controlled by pytest command-line option)
    if request.config.getoption("--cleanup"):
        print(f"Cleaning up {test_suite_dir}")
        shutil.rmtree(test_suite_dir)


def pytest_generate_tests(metafunc):
    if "test_suite_setup" in metafunc.fixturenames and "test_case" in metafunc.fixturenames:
        # Get the YAML file path from the command-line option
        yaml_file = metafunc.config.getoption("--test-config-file")
        test_suites = load_test_suites(yaml_file)

        # Parameterize tests
        params = []
        for suite in test_suites:
            if "tests" not in suite:
                pytest.fail(f"Test suite '{suite.get('name', 'unknown')}' missing 'tests' key.")
            for test in suite["tests"]:
                params.append((suite, test))
        
        metafunc.parametrize(
            ("test_suite_setup", "test_case"),
            params,
            indirect=["test_suite_setup"],
            ids=[f"{suite['name']}_{test['name']}" for suite, test in params],
        )

# Main test function
def test_converter(test_suite_setup, test_case):
    suite_name = test_suite_setup["suite_name"]
    converter_dir = test_suite_setup["converter_dir"]
    testdata_dir = test_suite_setup["testdata_dir"]
    testoutput_dir = test_suite_setup["testoutput_dir"]
    testconfig_dir = test_suite_setup["testconfig_dir"]
    testresults_dir = test_suite_setup["testresults_dir"]

    # print(f"suite_name = {suite_name}")
    # print(f"testdata_dir = {testdata_dir}")
    # print(f"testoutput_dir = {testoutput_dir}")
    # print(f"testconfig_dir = {testconfig_dir}")
    # print(f"testresults_dir = {testresults_dir}")
    # print(f"converter_dir = {converter_dir}")

    test_name = test_case["name"]
    converter = test_case["converter"]
    config_file = test_case.get("config")
    input_file = test_case.get("input")
    reference_file = test_case["reference"]

    # print(f"test_name = {test_name}")
    # print(f"converter = {converter}")
    # print(f"input_file = {input_file}")
    # print(f"reference_file = {reference_file}")
    # print(f"config_file = {config_file}")

    # Determine converter path
    if converter_dir:
        converter_path = Path(converter_dir) / converter
    else:
        converter_path = Path(converter)

    # Check if converter exists
    assert converter_path.exists(), f"Converter {converter_path} does not exist"

    # Check if required files exist
    if input_file:
        input_path = testdata_dir / input_file
        assert input_path.exists(), f"Input file {input_path} does not exist"
    if config_file:
        config_path = testconfig_dir / config_file
        assert config_path.exists(), f"Config file {config_path} does not exist"

    # output file name is assumed to be the same as the reference
    # file name
    output_file = reference_file
    output_path = testresults_dir / reference_file

    # Build and run the command
    cmd = [str(converter_path)]
    if config_file:
        cmd.extend(["--config", str(config_path)])
    if input_file and output_file:
        cmd.extend(["--input", str(input_path), "--output", str(output_path)])
    else:
        # If no input, assume config_file specifies input/output
        pytest.fail(f'Input/output or config specification required for {converter_path}')

    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    # Check if output file was created
    assert output_path.exists(), f"Output file {output_path} was not created"
    # could check reference earlier, but this allows the user to run
    # the converter and see that it generates output
    reference_path = testoutput_dir / reference_file
    assert reference_path.exists(), f"Reference file {reference_path} does not exist"

    # Compare output with reference
    assert run_nccmp(output_path, reference_path), f"Output {output_path} does not match reference {reference_path}"

if __name__ == "__main__":
    pytest.main(["-v", "--tb=short"])
