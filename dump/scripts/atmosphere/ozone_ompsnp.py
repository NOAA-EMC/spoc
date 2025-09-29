#!/usr/bin/env python3
import os
import numpy as np
import numpy.ma as ma

import bufr
from bufr.bufr_python.encoders import netcdf
from bufr.obs_builder import ObsBuilder, add_main_functions, add_dummy_variable, map_path
from bufr_ozone_obs_builder import BufrOzoneObsBuilder

MAPPING_PATH = map_path('bufr_ozone_ompsnp.yaml')
USE_REFERENCE_PRESSURE = True


class BufrOzoneOmpsnpObsBuilder(BufrOzoneObsBuilder):
    """
    Class for building observations from ompsn8 BUFR data.

    This class extends `ObsBuilder` to include specific logic for processing
    Level-2 retrived profile ozone data from OMPS nadir profiler

    :param mapping_path: Path to the mapping file.
    :type mapping_path: str
    """

    def __init__(self):
        super().__init__(MAPPING_PATH, log_name=os.path.basename(__file__))

    def make_obs(self, comm, input_path):

        # Get container from mapping file first
        self.log.info('Get container from bufr')
        container = super().make_obs(comm, input_path)

        self.log.debug(f'container list (original): {container.list()}')
        self.log.debug(f'all_sub_categories =  {container.all_sub_categories()}')
        self.log.debug(f'category map =  {container.get_category_map()}')

        # Add new/derived data into container
        for cat in container.all_sub_categories():

            self.log.debug(f'category = {cat}')

            satId = container.get('satelliteId', cat)
            if not np.any(satId):
                self.log.warning(f'category {cat[0]} does not exist in input file')

            self._add_ompsnp_pressures(container, cat)

        # Check
        self.log.debug(f'container list (updated): {container.list()}')
        self.log.debug(f'all_sub_categories {container.all_sub_categories()}')

        return container

    def _make_description(self):
        description = super()._make_description()
        self._add_new_pressure_descriptions(description)
        self._add_new_pressure_dimensions(description)

        return description

    def _add_new_pressure_dimensions(self, description):
        description = super()._make_description()
        description.add_dimension("Vertice", ["*/VERTICE"])

        return description

    def _add_new_pressure_descriptions(self, description):
        description.add_variables([
            {
                'name': 'RetrievalAncillaryData/pressureVertice',
                'source': 'pressureVertice',
                'units': 'Pa',
                'longName': 'Retrieval Pressure Vertices',
            }])

    def _add_ompsnp_pressures(self, container, cat):
        """
        Add pressure-related variables and data to the container.

        This includes generating pressure vertices and replacing the pressure variable
        as needed, depending on whether reference pressures are used.

        :param container: Observation container to modify.
        :type container: object
        :param cat: Category key for container access.
        :type cat: any
        """

        # Add new MetaData pressure-related variables
        satId = container.get('satelliteId', cat)
        if not satId.size:
            self.log.warning(f'category {cat[0]} does not exist in input file')
            pbot = container.get('bottomLevelPressure', cat)
            ptop = container.get('topLevelPressure', cat)
            pressureVertice = np.column_stack((ptop, pbot))
            pathLocation = container.get_paths('latitude', cat)
            pathPressure = [pathLocation[0], '*/VERTICE']
            container.add('pressureVertice', pressureVertice, pathPressure, cat)
            return

        o3val = container.get('ozoneLayer', cat)
        nlevs = 21
        nprofs = int(o3val.shape[0]/nlevs)

        # Set reference pressure
        if USE_REFERENCE_PRESSURE:
            # Set reference pressure vertices (define pressyre boundaries (in hPa))
            # These pressure vertices are defined in GSI ozinfo
            ptop = np.array([631.000, 398.000, 251.000, 158.000, 100.000, 63.100,
                             39.800, 25.100, 15.800, 10.000, 6.310, 3.980, 2.510,
                             1.580, 1.000, 0.631, 0.398, 0.251, 0.158, 0.100, 0])
            pbot = np.array([1000.000, 631.000, 398.000, 251.000, 158.000, 100.000,
                             63.100, 39.800, 25.100, 15.800, 10.000, 6.310, 3.980,
                             2.510, 1.580, 1.000, 0.631, 0.398, 0.251, 0.158, 0.100])

            # Convert to Pascals and account for standard atmospheric pressure conversion
            ptop = ptop * 100. * 1.01325
            pbot = pbot * 100. * 1.01325
            pref = pbot
            pressure = np.tile(pref, nprofs).astype(np.float32)

            # Create 2D pressure boundaries array
            pres = np.zeros((nlevs, 2), dtype=np.float32)
            pres[:, 0] = ptop
            pres[:, 1] = pbot
            pressureVertice = np.tile(pres, (nprofs, 1)).astype(np.float32)
        else:
            # Get pressure vertices
            pbot = container.get('bottomLevelPressure', cat)
            ptop = container.get('topLevelPressure', cat)

            # Assign reference pressures
            pref = pbot
            pressure = np.tile(pref, nprofs).astype(np.float32)

            # Create 2D pressure boundaries array
            pressureVertice = np.column_stack((ptop, pbot))

        self.log.debug(f'shape(pressureVertice) = {pressureVertice.shape}')
        self.log.debug(f'pressureVertice min/max = {pressureVertice.min()} {pressureVertice.max()}')

        # Add 2D pressure boundaries array to container
        pathLocation = container.get_paths('latitude', cat)
        pathPressure = [pathLocation[0], '*/VERTICE']
        container.add('pressureVertice', pressureVertice, pathPressure, cat)

        # Add 2D pressure to container
        pressure_original = container.get('pressure', cat)
        self.log.debug(f'pressure_original min/max = {pressure_original.min()} {pressure_original.max()}')
        self.log.debug(f'shape(pressure_original) = {pressure_original.shape}')
        self.log.debug(f'pressure min/max = {pressure.min()} {pressure.max()}')
        self.log.debug(f'shape(pressure) = {pressure.shape}')
        container.replace('pressure', pressure, cat)


# Add main functions create_obs_file or create_obs_group
add_main_functions(BufrOzoneOmpsnpObsBuilder)
