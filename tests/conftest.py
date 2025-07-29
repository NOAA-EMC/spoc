import os
import tarfile
import urllib.request
import io
import shutil
import pytest

DATA_URL = "https://ftp.emc.ncep.noaa.gov/static_files/public/spoc/spoc-0.0.0.tgz"

@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(os.path.join(request.fspath.dirname, 'data'))

@pytest.fixture(scope="session", autouse=True)
def download_test_data():
    """Download and extract required test data."""
    tests_dir = os.path.dirname(__file__)
    data_dir = os.path.join(tests_dir, "data")

    if os.path.isdir(data_dir) and os.listdir(data_dir):
        return data_dir

    downloaded_dir = os.path.join(tests_dir, "remote_data")

    try:
        with urllib.request.urlopen(DATA_URL) as response:
            archive_data = response.read()
        with tarfile.open(fileobj=io.BytesIO(archive_data), mode="r:gz") as tar:
            tar.extractall(path=tests_dir, filter="data")
    except Exception as exc:
        shutil.rmtree(downloaded_dir, ignore_errors=True)
        pytest.skip(f"Could not download test data: {exc}")

    data_dir = os.path.join(tests_dir, "data")
    shutil.move(downloaded_dir, data_dir)
    os.makedirs(os.path.join(data_dir, 'testrun'), exist_ok=True)

    return data_dir
