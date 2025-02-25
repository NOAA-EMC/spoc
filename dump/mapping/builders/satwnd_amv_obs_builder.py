import numpy as np

from bufr.obs_builder import ObsBuilder

class SatWndAmvObsBuilder(ObsBuilder):
    def __init__(self, input_path, mapping_path):
        super().__init__(input_path, mapping_path)

    def compute_wind_components(self, wdir, wspd):
        """
        Compute the U and V wind components from wind direction and wind speed.

        Parameters:
            wdir (array-like): Wind direction in degrees (meteorological convention: 0° = North, 90° = East).
            wspd (array-like): Wind speed.

        Returns:
            tuple: U and V wind components as numpy arrays with dtype float32.
        """
        wdir_rad = np.radians(wdir)  # Convert degrees to radians
        u = -wspd * np.sin(wdir_rad)
        v = -wspd * np.cos(wdir_rad)

        return u.astype(np.float32), v.astype(np.float32)
