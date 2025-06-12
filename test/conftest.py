import pytest

# Custom pytest option for cleanup
def pytest_addoption(parser):
    parser.addoption(
        "--cleanup",
        action="store_true",
        default=False,
        help="Clean up downloaded data and test results after tests",
    )

    parser.addoption(
        "--test-config-file",
        action="store",
        default="spoc-tests.yaml",
        help="Path to the YAML file containing test suites (default: spoc-tests.yaml.yaml)"
    )

