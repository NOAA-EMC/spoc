#!/usr/bin/env python3
import sys
import os
import argparse
import time
import numpy as np

from wxflow import Logger

import bufr
from bufr.obs_builder.logger import logging
from builders.satwnd_amv_obs_builder import SatWndAmvObsBuilder


log_level = os.getenv('LOG_LEVEL', 'INFO')
logger = Logger('BUFR2IODA_satwnd_amv_avhrr.py', level=log_level, colored_log=False)


class SatWndAmvAvhrrObsBuilder(SatWndAmvObsBuilder):
    def __init__(self, input_path, mapping_path):
        super().__init__(input_path, mapping_path)

    def _make_description(self):
        description = bufr.encoders.Description(self.mapping_path)

        description.add_variables([
            {
                'name': 'ObsType/windEastward',
                'source': 'obstype_uwind',
                'units': '1',
                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
            },
            {
                'name': 'ObsType/windNorthward',
                'source': 'obstype_vwind',
                'units': '1',
                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
            },
            {
                'name': 'ObsValue/windEastward',
                'source': 'windEastward',
                'units': 'm s-1',
                'longName': 'Eastward Wind Component',
            },
            {
                'name': 'ObsValue/windNorthward',
                'source': 'windNorthward',
                'units': 'm s-1',
                'longName': 'Northward Wind Component',
            },
            # MetaData/windGeneratingApplication will be inferred from generatingApplication
            # following a search for the proper generatingApplication column
            {
                'name': 'MetaData/windGeneratingApplication',
                'source': 'windGeneratingApplication',
                'units': '1',
                'longName': 'Wind Generating Application',
            },
            # MetaData/qualityInformationWithoutForecast will be inferred from qualityInformation
            # following a search for the proper generatingApplication column
            {
                'name': 'MetaData/qualityInformationWithoutForecast',
                'source': 'qualityInformationWithoutForecast',
                'units': 'percent',
                'longName': 'Quality Information Without Forecast',
            }
        ])

        return description

    def _get_QualityInformation_and_GeneratingApplication(self, comm, gnap2D, pccf2D, satID):
        # For METOP-A/B/C AVHRR data (satID 3,4,5), qi w/o forecast (qifn) is
        # packaged in same vector of qi with ga = 5 (QI without forecast), and EE
        # is packaged in same vector of qi with ga=7 (Estimated Error (EE) in m/s
        # converted to a percent confidence) shape (4,nobs).
        #
        # For NOAA-15/18/19 AVHRR data (satID 206,209,223), qi w/o forecast
        # (qifn) is packaged in same vector of qi with ga = 1 (EUMETSAT QI
        # without forecast), and EE is packaged in same vector of qi with ga=4
        # (Estimated Error (EE) in m/s converted to a percent confidence) shape
        # (4,nobs).
        #
        # Must conduct a search and extract the correct vector for gnap and qi
        # 0. Define the appropriate QI and EE search values, based on satID
        if np.all(np.isin(satID, [206, 209, 223])):  # NESDIS AVHRR set
            findQI = 1
            findEE = 4
        elif np.all(np.isin(satID, [3, 4, 5])):  # EUMETSAT AVHRR set
            findQI = 5
            findEE = 7
            # There is a catch: prior to 2023 AVHRR winds from EUMETSAT were formatted in the same
            # way as NESDIS AVHRR winds and both were passed through the NC005080 tank as a single
            # dataset. In that case, we need to actually set findQI=1 and findEE=4 here.
            # Let's do a preliminary check to see if any gnap2D values match findQI. If not, let's
            # automatically switch to findQI=1, findEE=4 and presume pre-2023 EUMETSAT AVHRR format
            if np.any(np.isin(gnap2D, [findQI])) == False:
                logging(comm, 'DEBUG',
                        f'NO GNAP VALUE OF {findQI} EXISTS FOR EUMETSAT AVHRR DATASET, PRESUMING PRE-2023 FORMATTING')
                findQI = 1
                findEE = 4
        else:
            logging(comm, 'DEBUG', f'satID set not found (all satID values follow):')
            for sid in np.unique(satID):
                logging(comm, 'DEBUG', f'satID: {sid}')
        logging(comm, 'DEBUG', f'BTH: findQI={findQI}')
        # 1. Find dimension-sizes of ga and qi (should be the same!)
        gDim1, gDim2 = np.shape(gnap2D)
        qDim1, qDim2 = np.shape(pccf2D)
        logging(comm, 'INFO', f'Generating Application and Quality Information SEARCH:')
        logging(comm, 'DEBUG', f'Dimension size of GNAP ({gDim1},{gDim2})')
        logging(comm, 'DEBUG', f'Dimension size of PCCF ({qDim1},{qDim2})')
        # 2. Initialize gnap and qifn as None, and search for dimension of
        #    ga with values of findQI. If the same column exists for qi, assign
        #    gnap to ga[:,i] and qifn to qi[:,i], else raise warning that no
        #    appropriate GNAP/PCCF combination was found
        gnap = None
        qifn = None
        for i in range(gDim2):
            if np.unique(gnap2D[:, i].squeeze()) == findQI:
                if i <= qDim2:
                    logging(comm, 'INFO', f'GNAP/PCCF found for column {i}')
                    gnap = gnap2D[:, i].squeeze()
                    qifn = pccf2D[:, i].squeeze()
                else:
                    logging(comm, 'INFO', f'ERROR: GNAP column {i} outside of PCCF dimension {qDim2}')
        if (gnap is None) & (qifn is None):
            raise ValueError(f'GNAP == {findQI} NOT FOUND OR OUT OF PCCF DIMENSION-RANGE, WILL FAIL!')
        # If EE is needed, key search on np.unique(gnap2D[:,i].squeeze()) == findEE instead
        # NOTE: Make sure to return np.float32 or np.int32 types as appropriate!!!
        return gnap.astype(np.int32), qifn.astype(np.int32)

    def _get_obs_type(self, swcm):
        """
        Determine the observation type based on `swcm` and `chanfreq`.

        Parameters:
            swcm (array-like): Switch mode values.
            chanfreq (array-like): Channel frequency values (Hz).

        Returns:
            numpy.ndarray: Observation type array.

        Raises:
            ValueError: If any `obstype` is unassigned.
        """

        obstype = swcm.copy()

        # Use numpy vectorized operations
        obstype = np.where(swcm == 1, 244, obstype)  # IRLW

        if not np.any(np.isin(obstype, [244])):
            raise ValueError("Error: Unassigned ObsType found ... ")

        return obstype.astype(np.int32)

    def _make_obs(self, comm):

        # Get container from mapping file first
        logging(comm, 'INFO', 'Get container from bufr')
        container = bufr.Parser(self.input_path, self.mapping_path).parse(comm)

        logging(comm, 'DEBUG', f'container list (original): {container.list()}')
        logging(comm, 'DEBUG', f'all_sub_categories =  {container.all_sub_categories()}')
        logging(comm, 'DEBUG', f'category map =  {container.get_category_map()}')

        # Add new/derived data into container
        for cat in container.all_sub_categories():

            logging(comm, 'DEBUG', f'category = {cat}')

            satid = container.get('satelliteId', cat)
            if satid.size == 0:
                logging(comm, 'WARNING', f'category {cat[0]} does not exist in input file')
                paths = container.get_paths('windComputationMethod', cat)
                obstype = container.get('windComputationMethod', cat)
                container.add('obstype_uwind', obstype, paths, cat)
                container.add('obstype_vwind', obstype, paths, cat)

                paths = container.get_paths('windSpeed', cat)
                wob = container.get('windSpeed', cat)
                container.add('windEastward', wob, paths, cat)
                container.add('windNorthward', wob, paths, cat)

                paths = container.get_paths('windComputationMethod', cat)
                dummy = container.get('windSpeed', cat)
                container.add('windGeneratingApplication', dummy, paths, cat)
                container.add('qualityInformationWithoutForecast', dummy, paths, cat)

            else:
                # Add new variables: ObsType/windEastward & ObsType/windNorthward
                swcm = container.get('windComputationMethod', cat)
                chanfreq = container.get('sensorCentralFrequency', cat)

                logging(comm, 'DEBUG', f'swcm min/max = {swcm.min()} {swcm.max()}')
                logging(comm, 'DEBUG', f'chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')

                obstype = self._get_obs_type(swcm)

                logging(comm, 'DEBUG', f'obstype = {obstype}')
                logging(comm, 'DEBUG', f'obstype min/max =  {obstype.min()} {obstype.max()}')

                paths = container.get_paths('windComputationMethod', cat)
                container.add('obstype_uwind', obstype, paths, cat)
                container.add('obstype_vwind', obstype, paths, cat)

                # Add new variables: ObsValue/windEastward & ObsValue/windNorthward
                wdir = container.get('windDirection', cat)
                wspd = container.get('windSpeed', cat)

                logging(comm, 'DEBUG', f'wdir min/max = {wdir.min()} {wdir.max()}')
                logging(comm, 'DEBUG', f'wspd min/max = {wspd.min()} {wspd.max()}')

                uob, vob = self.compute_wind_components(wdir, wspd)

                logging(comm, 'DEBUG', f'uob min/max = {uob.min()} {uob.max()}')
                logging(comm, 'DEBUG', f'vob min/max = {vob.min()} {vob.max()}')

                paths = container.get_paths('windSpeed', cat)
                container.add('windEastward', uob, paths, cat)
                container.add('windNorthward', vob, paths, cat)

                # Add new variables: MetaData/windGeneratingApplication and qualityInformationWithoutForecast
                satID = container.get('satelliteId', cat)
                gnap2D = container.get('generatingApplication', cat)
                pccf2D = container.get('qualityInformation', cat)

                gnap, qifn = self._get_QualityInformation_and_GeneratingApplication(comm, gnap2D, pccf2D, satID)

                logging(comm, 'DEBUG', f'gnap min/max = {gnap.min()} {gnap.max()}')
                logging(comm, 'DEBUG', f'qifn min/max = {qifn.min()} {qifn.max()}')

                paths = container.get_paths('windComputationMethod', cat)
                container.add('windGeneratingApplication', gnap, paths, cat)
                container.add('qualityInformationWithoutForecast', qifn, paths, cat)

        # Check
        logging(comm, 'DEBUG', f'container list (updated): {container.list()}')
        logging(comm, 'DEBUG', f'all_sub_categories {container.all_sub_categories()}')

        return container


def create_obs_group(input_path, mapping_path, category, env):
    obs_builder = SatWndAmvAvhrrObsBuilder(input_path, mapping_path)
    obs_builder.create_obs_group(category, env)


def create_obs_file(input_path, mapping_path, output_path, type='netcdf', append=False):
    obs_builder = SatWndAmvAvhrrObsBuilder(input_path, mapping_path)
    return obs_builder.create_obs_file(output_path, type, append)


if __name__ == '__main__':
    start_time = time.time()

    bufr.mpi.App(sys.argv)
    comm = bufr.mpi.Comm("world")

    # Required input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input BUFR')
    parser.add_argument('mapping', type=str, help='BUFR2IODA Mapping File')
    parser.add_argument('output', type=str, help='Output NetCDF')

    args = parser.parse_args()
    create_obs_file(args.input, args.mapping, args.output)

    end_time = time.time()
    running_time = end_time - start_time
    logging(comm, 'INFO', f'Total running time: {running_time}')
