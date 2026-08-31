import sys, os
from pathlib import Path

def config_path(*path_components):
    config_path = os.path.realpath(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..', 'config', *path_components
        )
    )

    if os.path.exists(config_path):
        return config_path
    else:
        raise FileNotFoundError(f"Configuration file not found: {os.path.join(*path_components)}")
