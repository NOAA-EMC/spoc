#!/usr/bin/env python3

import os
import numpy as np
import numpy.ma as ma
from pathlib import Path
from typing import Dict


import bufr
from bufr.obs_builder import ObsBuilder

class BaseGpsroBufrObsBuilder(ObsBuilder):
    """Common logic; subclasses supply the YAML map paths."""

    def __init__(self, map_dict: Dict[str, str], *, log_name: str):
        # keep a copy so subclasses can use it later
        self.map_dict: Dict[str, str] = {k: str(v) for k, v in map_dict.items()}

        # promote the two paths to attributes for easy access
        self.loc_profile: str = self.map_dict["loc_profile"]
        self.height_profile: str = self.map_dict["height_profile"]

        # fast sanity-check
        for p in (self.loc_profile, self.height_profile):
            if not Path(p).is_file():
                raise FileNotFoundError(f"Mapping file not found: {p}")

        # call the real ObsBuilder ctor
        super().__init__(self.map_dict, log_name=log_name)

    # -----------------------------------------------------------------
    # Example method the parent framework calls
    # -----------------------------------------------------------------
    def make_obs(self,  comm, input_path):
        """
        Create the ioda gpsro bufr observations:
        - reads values
        - adds sequenceNum

        Parameters
        ----------
        comm: object
                The communicator object (e.g., MPI)
        input_path: str
                The input bufr file
        """

        self.log.info("Get the two separate latitude and height containers")
        container = bufr.Parser(input_path, self.loc_profile).parse(comm)
        height_container = bufr.Parser(input_path, self.height_profile).parse(comm)

        # The height and location profiles are the same shape, with each height profile corresponding
        # to the location profile. So lets merge them together.
        self.log.info("Merge the separate location and height containers")
        for cat in container.all_sub_categories():
            for var in height_container.list():
                container.add(var, height_container.get(var, cat), ['*'], cat)

        self.log.debug(f'container list (original): {container.list()}') #['atmosphericRefractivity', 'bendingAngle_roseq2repl1']
        self.log.debug(f'all_sub_categories =  {container.all_sub_categories()}') # [['cosmic2_750'], ['cosmic2_751']]
        self.log.debug(f'category map =  {container.get_category_map()}') # {'splits/satId': ['metop_3', 'metop_4']}

        self.log.debug(f'Create new datacontainer')
        d1 = bufr.DataContainer({'splits/satId': ['metop', 'cosmic', 'tdm', 'grace', 'geoopt', 'piq', 'k5', 'paz', 's6', 'tsx', 'spire']})
        #d1 = bufr.DataContainer({'splits/satId': ['metop', 'cosmic', 'tdm']}) #, 'grace', 'geoopt', 'piq', 'k5', 'paz', 's6', 'tsx', 'spire']})

        # Lists required.
        metop_list = [['metop_3'], ['metop_4'], ['metop_5']]
        metop_variables = {}
        metop_paths = {}
        cosmic_list = [['cosmic2_750'], ['cosmic2_751'], ['cosmic2_752'], ['cosmic2_753'], ['cosmic2_754'], ['cosmic2_755']]
        cosmic_variables = {}
        cosmic_paths = {}
        tdm_list = [['tdm_43']]
        tdm_variables = {}
        tdm_paths = {}
        grace_list = [['grace_803'], ['grace_804']]
        grace_variables = {}
        grace_paths = {}
        geoopt_list = [['geoopt_265'], ['geoopt_266']]
        geoopt_variables = {}
        geoopt_paths = {}
        piq_list = [['piq_267'], ['piq_268']]
        piq_variables = {}
        piq_paths = {}
        k5_list = [['k5_825']]
        k5_variables = {}
        k5_paths = {}
        paz_list = [['paz_44']]
        paz_variables = {}
        paz_paths = {}
        s6_list = [['s6_66']]
        s6_variables = {}
        s6_paths = {}
        tdm_list = [['tdm_43']]
        tdm_variables = {}
        tdm_paths = {}
        tsx_list = [['tsx_42']]
        tsx_variables = {}
        tsx_paths = {}
        spire_list = [['spire_269']]
        spire_variables = {}
        spire_paths = {}

        self.log.info("Process by category")
        number_satellites_processed = 0
        for cat in container.all_sub_categories():

            self.log.debug(f'category = {cat}')

            satellite_variables = {}
            satellite_paths = {}

            satId = container.get('satelliteId', cat)
            if not np.any(satId):
                self.log.warning(f'category {cat[0]} does not exist in input file')
            else:
                self.log.info(f"For {cat}, retrieve data from container")
                if len(container.list()) == 0:
                    self.log.warning(f'category {cat} does not contain any data')
                else:

                    # do manipulations first and add them to the container
                    self.log.info("   Creating derived variables - stationIdentification")
                    self._derive_stationidentification(container, cat)

                    self.log.info("   Creating and replacing derived variables - Grid Latitude / Longitude")
                    self._replace_gridcoordinates(container, cat)

                    self.log.info("   Creating derived variables - impactHeightRO")
                    self._derive_imph(container, cat)

                    self.log.info("   Updating bnda, imph, impp, bndaoe depending on mefr")
                    self._update_depending_on_mefr(container, cat)

                    self.log.info("   Updating Sequence Number")
                    self._update_sequencenumber(container, cat)

                    self.log.info("   Updating satelliteAscendingFlag and QFRO")
                    self._update_satelliteascendingflag_and_qualityflags(container, cat)


                    # GLOBAL ATTRIBUTES


                    # Get generic information for arrays for empty output files
                    self.log.info("   - Gather information for empty output files")
                    if number_satellites_processed < 1:
                        varlist = np.array(container.list())
                        dtypes = {}
                        for varnames in varlist:
                            dtypes[varnames] = container.get(varnames, cat).dtype
                            #print(f"NICKE dtypes {varnames}: {dtypes[varnames]}")
                        number_satellites_processed+=1

                    self.log.info("   - Get data from container and separate by group (metop, cosmic, etc)")
                    # THEN look in the container, put everything in arrays, and do everything.
                    for varname in container.list():  # varname is a string
                        satellite_variables[varname]= np.array(container.get(varname, cat))
                        satellite_paths[varname] = container.get_paths(varname, cat)

                        if cat in metop_list:
                            if varname not in metop_variables:
                                metop_variables[varname] = satellite_variables[varname]
                                metop_paths[varname] = satellite_paths[varname]
                            else:
                                metop_variables[varname] = np.concatenate((satellite_variables[varname], metop_variables[varname]), axis=0)
                                metop_paths[varname] = satellite_paths[varname]

                        elif cat in cosmic_list:
                            if varname not in cosmic_variables:
                                cosmic_variables[varname] = satellite_variables[varname]
                                cosmic_paths[varname] = satellite_paths[varname]
                            else:
                                cosmic_variables[varname] = np.concatenate((satellite_variables[varname], cosmic_variables[varname]), axis=0)
                                cosmic_paths[varname] = satellite_paths[varname]

                        elif cat in tdm_list:
                            if varname not in tdm_variables:
                                tdm_variables[varname] = satellite_variables[varname]
                                tdm_paths[varname] = satellite_paths[varname]
                            else:
                                tdm_variables[varname] = np.concatenate((satellite_variables[varname], tdm_variables[varname]), axis=0)
                                tdm_paths[varname] = satellite_paths[varname]

                        elif cat in grace_list:
                            if varname not in grace_variables:
                                grace_variables[varname] = satellite_variables[varname]
                                grace_paths[varname] = satellite_paths[varname]
                            else:
                                grace_variables[varname] = np.concatenate((satellite_variables[varname], grace_variables[varname]), axis=0)
                                grace_paths[varname] = satellite_paths[varname]

                        elif cat in geoopt_list:
                            if varname not in geoopt_variables:
                                geoopt_variables[varname] = satellite_variables[varname]
                                geoopt_paths[varname] = satellite_paths[varname]
                            else:
                                geoopt_variables[varname] = np.concatenate((satellite_variables[varname], geoopt_variables[varname]), axis=0)
                                geoopt_paths[varname] = satellite_paths[varname]

                        elif cat in piq_list:
                            if varname not in piq_variables:
                                piq_variables[varname] = satellite_variables[varname]
                                piq_paths[varname] = satellite_paths[varname]
                            else:
                                piq_variables[varname] = np.concatenate((satellite_variables[varname], piq_variables[varname]), axis=0)
                                piq_paths[varname] = satellite_paths[varname]

                        elif cat in k5_list:
                            if varname not in k5_variables:
                                k5_variables[varname] = satellite_variables[varname]
                                k5_paths[varname] = satellite_paths[varname]
                            else:
                                k5_variables[varname] = np.concatenate((satellite_variables[varname], k5_variables[varname]), axis=0)
                                k5_paths[varname] = satellite_paths[varname]

                        elif cat in paz_list:
                            if varname not in paz_variables:
                                paz_variables[varname] = satellite_variables[varname]
                                paz_paths[varname] = satellite_paths[varname]
                            else:
                                paz_variables[varname] = np.concatenate((satellite_variables[varname], paz_variables[varname]), axis=0)
                                paz_paths[varname] = satellite_paths[varname]

                        elif cat in s6_list:
                            if varname not in s6_variables:
                                s6_variables[varname] = satellite_variables[varname]
                                s6_paths[varname] = satellite_paths[varname]
                            else:
                                s6_variables[varname] = np.concatenate((satellite_variables[varname], s6_variables[varname]), axis=0)
                                s6_paths[varname] = satellite_paths[varname]

                        elif cat in tsx_list:
                            if varname not in tsx_variables:
                                tsx_variables[varname] = satellite_variables[varname]
                                tsx_paths[varname] = satellite_paths[varname]
                            else:
                                tsx_variables[varname] = np.concatenate((satellite_variables[varname], tsx_variables[varname]), axis=0)
                                tsx_paths[varname] = satellite_paths[varname]

                        elif cat in spire_list:
                            if varname not in spire_variables:
                                spire_variables[varname] = satellite_variables[varname]
                                spire_paths[varname] = satellite_paths[varname]
                            else:
                                spire_variables[varname] = np.concatenate((satellite_variables[varname], spire_variables[varname]), axis=0)
                                spire_paths[varname] = satellite_paths[varname]



        # Needs to stay this way.
        # Once manipulations are done and all the concantenations are done, THEN you can put all the data where they need to be.
        self.log.info("Add all data to new container.")
        if metop_variables:
            self.log.info(" - Adding metop")
            for varname, data in metop_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), metop_paths[varname], ['metop'])
        else:
            self.log.info(" - No Metop variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['metop'])

        if cosmic_variables:
            self.log.info(" - Adding cosmic")
            for varname, data in cosmic_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), cosmic_paths[varname], ['cosmic'])
        else:
            self.log.info(" - No Cosmic variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['cosmic'])

        if tdm_variables:
            self.log.info(" - Adding tdm")
            for varname, data in tdm_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), tdm_paths[varname], ['tdm'])
        else:
            self.log.info(" - No tdm variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['tdm'])

        if grace_variables:
            for varname, data in grace_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), grace_paths[varname], ['grace'])
        else:
            self.log.info("No grace variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['grace'])

        if geoopt_variables:
            for varname, data in geoopt_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), geoopt_paths[varname], ['geoopt'])
        else:
            self.log.info("No geoopt variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['geoopt'])

        if piq_variables:
            for varname, data in piq_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), piq_paths[varname], ['piq'])
        else:
            self.log.info("No piq variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['piq'])

        if k5_variables:
            for varname, data in k5_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), k5_paths[varname], ['k5'])
        else:
            self.log.info("No k5 variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['k5'])

        if paz_variables:
            for varname, data in paz_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), paz_paths[varname], ['paz'])
        else:
            self.log.info("No paz variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['paz'])

        if s6_variables:
            for varname, data in s6_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), s6_paths[varname], ['s6'])
        else:
            self.log.info("No s6 variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['s6'])

        if tsx_variables:
            for varname, data in tsx_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), tsx_paths[varname], ['tsx'])
        else:
            self.log.info("No tsx variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['tsx'])

        if spire_variables:
            for varname, data in spire_variables.items():
                self.log.debug(f"      Variable and shape: {varname} {data.shape}")
                d1.add(str(varname), np.array(data), spire_paths[varname], ['spire'])
        else:
            self.log.info("No spire variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['spire'])

        # Check
        self.log.debug(f'container list (updated): {container.list()}')

        #return container
        return d1


    # Provide defualt implementations for methods from the ObsBuilder class
    def _make_description(self):
        description = super()._make_description()
        self._add_new_variable_descriptions(description)
        self._add_new_variables_height_yaml(description)

        return description

    def _add_new_variable_descriptions(self, description):
        description.add_variables([
            {
                'name': 'MetaData/stationIdentification',
                'source': 'stationIdentification',
                'longName': 'Station Identification',
            },
            {
                'name': 'MetaData/impactHeightRO',
                'source': 'impactHeightRO1',
                'units': 'm',
                'longName': 'Impact Height Bending Angle',
            }
            ])

    def _add_new_variables_height_yaml(self, description):
        description.add_variables([
            {
                'name': 'MetaData/height',
                'source': 'height',
                'units': 'm',
                'longName': 'Height for Atm Refractivity',
            },
            {
                'name': 'MetaData/percentConfidence',
                'source': 'percentConfidence',
                'units': '%',
                'longName': 'Ref Percent Confidence',
            },
            {
                'name': 'ObsValue/atmosphericRefractivity',
                'source': 'atmosphericRefractivity',
                'units': 'N-units',
                'longName': 'Atmospheric Refractivity',
            },
            {
                'name': 'ObsError/atmosphericRefractivity',
                'source': 'obsErrorAtmosphericRefractivity',
                'units': 'N-units',
                'longName': 'Atmospheric Refractivity Obs Error',
            },
            ])


    # Methods that are used to extend the export description
    def _replace_gridcoordinates(self, container, cat):
        lat_deg = container.get('latitude', cat)
        lon_deg = container.get('longitude', cat)

        lat_rad = np.deg2rad(lat_deg)
        lon_rad = np.deg2rad(lon_deg)

        # Add to container
        container.replace('gridLatitude', lat_rad, cat)
        container.replace('gridLongitude', lon_rad, cat)


    def _derive_stationidentification(self, container, cat):
        said = container.get('satelliteId', cat)
        ptid = container.get('satelliteTransmitterId', cat)
        said_paths = container.get_paths('satelliteId', cat)

        stid = []
        for i in range(len(said)):
            newval = str(said[i]).zfill(4)+str(ptid[i]).zfill(4)
            stid.append(str(newval))
        stid = np.array(stid).astype(dtype='str')
        stid = ma.array(stid)
        ma.set_fill_value(stid, "")

        # Add to container
        container.add('stationIdentification', stid, said_paths, cat)


    def _derive_imph(self, container, cat):
        impp1 = container.get('impactParameterRO_roseq2repl1', cat).astype(np.float32)
        impp2 = container.get('impactParameterRO_roseq2repl2', cat).astype(np.float32)
        impp3 = container.get('impactParameterRO_roseq2repl3', cat).astype(np.float32)
        elrc = container.get('earthRadiusCurvature', cat)
        geodu = container.get('geoidUndulation', cat)
        impp1_paths = container.get_paths('impactParameterRO_roseq2repl1', cat)

        # Calculate imph
        imph1 = (impp1 - elrc - geodu).astype(np.float32)
        imph2 = (impp2 - elrc - geodu).astype(np.float32)
        imph3 = (impp3 - elrc - geodu).astype(np.float32)

        # Add to container
        container.add('impactHeightRO1', imph1, impp1_paths, cat)
        container.add('impactHeightRO2', imph2, impp1_paths, cat)
        container.add('impactHeightRO3', imph3, impp1_paths, cat)


    def _update_depending_on_mefr(self, container, cat):
        mefr1 = container.get('frequency__roseq2repl1', cat)
        mefr2 = container.get('frequency__roseq2repl2', cat)
        mefr3 = container.get('frequency__roseq2repl3', cat)
        impp1 = container.get('impactParameterRO_roseq2repl1', cat).astype(np.float32)
        impp2 = container.get('impactParameterRO_roseq2repl2', cat).astype(np.float32)
        impp3 = container.get('impactParameterRO_roseq2repl3', cat).astype(np.float32)
        bnda1 = container.get('bendingAngle_roseq2repl1', cat)
        bnda2 = container.get('bendingAngle_roseq2repl2', cat)
        bnda3 = container.get('bendingAngle_roseq2repl3', cat)
        bndaoe1 = container.get('obsErrorBendingAngle1', cat)
        bndaoe2 = container.get('obsErrorBendingAngle2', cat)
        bndaoe3 = container.get('obsErrorBendingAngle3', cat)
        imph1 = container.get('impactHeightRO1', cat)
        imph2 = container.get('impactHeightRO2', cat)
        imph3 = container.get('impactHeightRO3', cat)

        self.log.info("      Overwrite values for mefr, bnda, impp, imph, bndaoe")
        for i in range(len(impp1)):
            if (mefr2[i] == 0.0):
                mefr1[i] = mefr2[i]
                bnda1[i] = bnda2[i]
                impp1[i] = impp2[i]
                imph1[i] = imph2[i]
                bndaoe1[i] = bndaoe2[i]
            if (mefr3[i] == 0.0):
                mefr1[i] = mefr3[i]
                bnda1[i] = bnda3[i]
                impp1[i] = impp3[i]
                imph1[i] = imph3[i]
                bndaoe1[i] = bndaoe3[i]

        # Replace in container
        container.replace('bendingAngle_roseq2repl1', bnda1, cat)
        container.replace('frequency__roseq2repl1', mefr1, cat)
        container.replace('impactParameterRO_roseq2repl1', impp1, cat)
        container.replace('impactHeightRO1', imph1, cat)
        container.replace('obsErrorBendingAngle1', bndaoe1, cat)


    def _update_sequencenumber(self, container, cat):
        seqnum = container.get('sequenceNumber', cat)

        count1 = 0
        count2 = 0
        seqnum2 = []
        for i in range(len(seqnum)):
            if (int(seqnum[i]) != count2):
                count1 += 1
            count2 = int(seqnum[i])
            seqnum2.append(count1)
        seqnum2 = np.array(seqnum2)

        # Add to container
        container.replace('sequenceNumber', seqnum2, cat)

    def _update_satelliteascendingflag_and_qualityflags(self, container, cat):
        qfro = container.get('qualityFlags', cat)
        qfro2 = container.get('pccf', cat).astype(np.float32)
        satasc = container.get('satelliteAscendingFlag', cat)
        #   find ibit for qfro (16bit from left to right)
        #   bit5=1, reject the bending angle obs
        #   bit6=1, reject the refractivity obs
        bit3 = []
        bit5 = []
        bit6 = []
        for quality in qfro:
            if quality & 8192 > 0:
                bit3.append(1)
            else:
                bit3.append(0)

            if quality & 2048 > 0:
                bit5.append(1)
            else:
                bit5.append(0)

            # For refractivity data use only:
            if quality & 1024 > 0:
                bit6.append(1)
            else:
                bit6.append(0)

        bit3 = np.array(bit3)
        bit5 = np.array(bit5)
        bit6 = np.array(bit6)

        # overwrite satelliteAscendingFlag and QFRO
        for quality in range(len(bit3)):
            satasc[quality] = 0
            qfro2[quality] = 0.0
            if bit3[quality] == 1:
                satasc[quality] = 1
            # if (bit6[quality] == 1): refractivity data only
            #    qfro2[quality] = 1.0
            if (bit5[quality] == 1):
                qfro2[quality] = 1.0

        # Add to container
        container.replace('qualityFlags', qfro2, cat)
        container.replace('satelliteAscendingFlag', satasc, cat)
