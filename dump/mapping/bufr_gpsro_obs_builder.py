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

        container = bufr.Parser(input_path, self.loc_profile).parse(comm)
        height_container = bufr.Parser(input_path, self.height_profile).parse(comm)
       
        # The height and location profiles are the same shape, with each height profile corresponding
        # to the location profile. So lets merge them together.
        self.log.info("Merge the separate location and height containers")
        for cat in container.all_sub_categories():
            for var in height_container.list():
                container.add(var, height_container.get(var, cat), ['*'], cat)
           

        self.log.info("Add variables")

        # Get container from mapping file first
        #container = bufr.Parser(input_path, MAPPING_PATH).parse(comm)

        self.log.debug(f'container list (original): {container.list()}') #['atmosphericRefractivity', 'bendingAngle_roseq2repl1']
        self.log.debug(f'all_sub_categories =  {container.all_sub_categories()}') # [['cosmic2_750'], ['cosmic2_751']]
        self.log.debug(f'category map =  {container.get_category_map()}') # {'splits/satId': ['metop_3', 'metop_4']}

##### NICKE get general list
        #varlist = np.array(container.list())
        #print("NICKE varlist ", varlist)
        #dtypes = {}
        #for varnames in varlist:
        #    dtypes[varnames] = container.get(varnames, ['cosmic2_750']).dtype
        #    print(f"NICKE dtypes {varnames}: {dtypes[varnames]}")

###### NICKE Make the new container
        self.log.debug(f'NNN create satId datacontainer')
        #d1 = bufr.DataContainer({'splits/satId': ['metop', 'cosmic', 'tdm', 'grace', 'geoopt', 'piq', 'k5', 'paz', 's6', 'tsx', 'spire']})
        d1 = bufr.DataContainer({'splits/satId': ['metop', 'cosmic', 'tdm']}) #, 'grace', 'geoopt', 'piq', 'k5', 'paz', 's6', 'tsx', 'spire']})

###### NICKE make lists
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

        #maybe figure out how to use this. idk how right now.
        #full_list = ['metop_list', 'cosmic_list', 'tdm_list']

        # Add new/derived data into container
        number_satellites_processed = 0
        for cat in container.all_sub_categories():

            self.log.debug(f'category = {cat}')

            satellite_variables = {}
            satellite_paths = {}

            satId = container.get('satelliteId', cat)
            if not np.any(satId):
                self.log.warning(f'category {cat[0]} does not exist in input file')
            else:     # if it exists, move it to a category 
##### NICKE BEGIN MERGE
                self.log.info(f"For {cat}, retrieve data from container")
                #number_satellites_processed = 0

                #get the data for the category
                if len(container.list()) == 0:
                    self.log.warning(f'category {cat} does not contain any data')
                else:

                    # do manipulations first and add them to the container
                    self.log.info("Manipulate the data")

                    self.log.info("Creating derived variables - stationIdentification")
                    self._derive_stationIdentification(container, cat)
                    #print("NICKE NEW CONTAINER 2 LIST", container.list())

                    self.log.info("Creating derived variables - Grid Latitude / Longitude")
                    #print(f"NICKE gridlat before {satellite_variables['gridLatitude'].min()}, {satellite_variables['gridLatitude'].max()}")
                    #self._replace_gridCoordinates(container, cat)   # Don't need this but I know it works if I get rid of conversion in yaml!


                    self.log.info("Deriving imph and using mefr to manipulate bnda, imph, impp, bndaoe")
                    #self._derive_imph_and_bnda(container, cat)


                    # Get generic information for arrays for empty output files 
                    self.log.info("Gather information for empty output files")
                    if number_satellites_processed < 1:
                        varlist = np.array(container.list())
                        #print("NICKE varlist ", varlist)
                        dtypes = {}
                        for varnames in varlist:
                            dtypes[varnames] = container.get(varnames, cat).dtype
                            #print(f"NICKE dtypes {varnames}: {dtypes[varnames]}")
                        number_satellites_processed+=1
                     

                    self.log.info("Get data from container and separate by group (metop, cosmic, etc)")
                    # THEN look in the container, put everything in arrays, and do everything.
                    for varname in container.list():  # varname is a string
                        #print(f"NICKE varname {varname}")
                        satellite_variables[varname]= np.array(container.get(varname, cat))
                        satellite_paths[varname] = container.get_paths(varname, cat)
                        #print(f"NICKE thedata {thedata.shape}, {thedata}")
                        #print(f"NICKE satellite_variables[varname]: {satellite_variables[varname]}")
   
    
                        if cat in metop_list: 
                            print(f"NICKE cat {cat}")
                            if varname not in metop_variables:
                                metop_variables[varname] = satellite_variables[varname] 
                                metop_paths[varname] = satellite_paths[varname] 
                                #print(f"NICKE metop ADD {cat} {varname}")
                            else:
                                metop_variables[varname] = np.concatenate((satellite_variables[varname], metop_variables[varname]), axis=0)
                                metop_paths[varname] = satellite_paths[varname]
                                #print(f"NICKE metop APPEND {cat} {varname}")
    
                        elif cat in cosmic_list: 
                            print(f"NICKE cat {cat}")
    
                            if varname not in cosmic_variables:
                                cosmic_variables[varname] = satellite_variables[varname]  # start a list of arrays
                                cosmic_paths[varname] = satellite_paths[varname]
                                #print(f"NICKE cosmic ADD {cat} {varname} {satellite_paths[varname]}")
                            else:
                                cosmic_variables[varname] = np.concatenate((satellite_variables[varname], cosmic_variables[varname]), axis=0)
                                cosmic_paths[varname] = satellite_paths[varname]
                                #print(f"NICKE cosmic APPEND {cat} {varname} {satellite_paths[varname]}")
    
    
                        elif cat in tdm_list:
                            print(f"NICKE cat {cat}")
    
                            if varname not in tdm_variables:
                                tdm_variables[varname] = satellite_variables[varname]  # start a list of arrays
                                tdm_paths[varname] = satellite_paths[varname]
                                #print(f"NICKE tdm ADD {cat} {varname} {satellite_paths[varname]}")
                            else:
                                tdm_variables[varname] = np.concatenate((satellite_variables[varname], tdm_variables[varname]), axis=0)
                                tdm_paths[varname] = satellite_paths[varname]
                                #print(f"NICKE tdm APPEND {cat} {varname} {satellite_paths[varname]}")
    

#                    elif cat in grace_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in grace_variables:
#                            grace_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            grace_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE grace ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            grace_variables[varname] = np.concatenate((satellite_variables[varname], grace_variables[varname]), axis=0)
#                            grace_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE grace APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in geoopt_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in geoopt_variables:
#                            geoopt_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            geoopt_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE geoopt ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            geoopt_variables[varname] = np.concatenate((satellite_variables[varname], geoopt_variables[varname]), axis=0)
#                            geoopt_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE geoopt APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in piq_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in piq_variables:
#                            piq_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            piq_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE piq ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            piq_variables[varname] = np.concatenate((satellite_variables[varname], piq_variables[varname]), axis=0)
#                            piq_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE piq APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in k5_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in k5_variables:
#                            k5_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            k5_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE k5 ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            k5_variables[varname] = np.concatenate((satellite_variables[varname], k5_variables[varname]), axis=0)
#                            k5_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE k5 APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in paz_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in paz_variables:
#                            paz_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            paz_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE paz ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            paz_variables[varname] = np.concatenate((satellite_variables[varname], paz_variables[varname]), axis=0)
#                            paz_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE paz APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in s6_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in s6_variables:
#                            s6_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            s6_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE s6 ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            s6_variables[varname] = np.concatenate((satellite_variables[varname], s6_variables[varname]), axis=0)
#                            s6_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE s6 APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in tsx_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in tsx_variables:
#                            tsx_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            tsx_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE tsx ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            tsx_variables[varname] = np.concatenate((satellite_variables[varname], tsx_variables[varname]), axis=0)
#                            tsx_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE tsx APPEND {cat} {varname} {satellite_paths[varname]}")
#
#
#                    elif cat in spire_list:
#                        print(f"NICKE cat {cat}")
#
#                        if varname not in spire_variables:
#                            spire_variables[varname] = satellite_variables[varname]  # start a list of arrays
#                            spire_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE spire ADD {cat} {varname} {satellite_paths[varname]}")
#                        else:
#                            spire_variables[varname] = np.concatenate((satellite_variables[varname], spire_variables[varname]), axis=0)
#                            spire_paths[varname] = satellite_paths[varname]
#                            print(f"NICKE spire APPEND {cat} {varname} {satellite_paths[varname]}")





        # Needs to stay this way.
        # Once manipulations are done and all the concantenations are done, THEN you can put all the data where they need to be.
        self.log.info("Add all data to new container.")
        if metop_variables:
            self.log.info(" - Adding metop")
            for varname, data in metop_variables.items():
                print(f"\nVariable: {varname}")
                print(f"Shape: {data.shape}")
                d1.add(str(varname), np.array(data), metop_paths[varname], ['metop'])
        else:
            print(" - No Metop variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['metop'])

        if cosmic_variables:
            self.log.info(" - Adding cosmic")
            for varname, data in cosmic_variables.items():
                #print(f"\nVariable: {varname}")
                #print(f"Shape: {data.shape}")
                d1.add(str(varname), np.array(data), cosmic_paths[varname], ['cosmic'])
        else:
            print(" - No Cosmic variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['cosmic'])

        if tdm_variables:
            self.log.info(" - Adding tdm")
            for varname, data in tdm_variables.items():
                #print(f"\nVariable: {varname}")
                #print(f"Shape: {data.shape}")
                d1.add(str(varname), np.array(data), tdm_paths[varname], ['tdm'])
        else:
            print(" - No tdm variables were collected. Output file will be empty.")
            for varname in varlist:
                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['tdm'])

#        if grace_variables:
#            for varname, data in grace_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), grace_paths[varname], ['grace'])
#        else:
#            print("No grace variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['grace'])
#           
#        if geoopt_variables:
#            for varname, data in geoopt_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), geoopt_paths[varname], ['geoopt'])
#        else:
#            print("No geoopt variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['geoopt'])
#
#        if piq_variables:
#            for varname, data in piq_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), piq_paths[varname], ['piq'])
#        else:
#            print("No piq variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['piq'])
#
#        if k5_variables:
#            for varname, data in k5_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), k5_paths[varname], ['k5'])
#        else:
#            print("No k5 variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['k5'])
#
#        if paz_variables:
#            for varname, data in paz_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), paz_paths[varname], ['paz'])
#        else:
#            print("No paz variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['paz'])
#
#        if s6_variables:
#            for varname, data in s6_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), s6_paths[varname], ['s6'])
#        else:
#            print("No s6 variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['s6'])
#
#        if tsx_variables:
#            for varname, data in tsx_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), tsx_paths[varname], ['tsx'])
#        else:
#            print("No tsx variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['tsx'])
#
#        if spire_variables:
#            for varname, data in spire_variables.items():
#                print(f"\nVariable: {varname}")
#                print(f"Shape: {data.shape}")
#                d1.add(str(varname), np.array(data), spire_paths[varname], ['spire'])
#        else:
#            print("No spire variables were collected. Output file will be empty.")
#            for varname in varlist:
#                d1.add(str(varname), np.array([]).astype(dtypes[varname]), ['*'], ['spire'])



##### NICKE END MERGE
        

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
                #'units': 'Pa',
                'longName': 'Station Identification',
            }])

    def _add_new_variables_height_yaml(self, description):
        description.add_variables([
            {
                'name': 'MetaData/satelliteId_height',
                'source': 'satelliteId_height',
                #'units': 'Pa',
                'longName': 'satelliteId_height',
            },
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
    #def _add_pressures(self, container, cat):
    def _replace_gridCoordinates(self, container, cat):
        lat_deg = container.get('gridLatitude', cat) 
        lon_deg = container.get('gridLongitude', cat) 
        lat_deg_paths = container.get_paths('gridLatitude', cat) 
        lon_deg_paths = container.get_paths('gridLongitude', cat)

        print(f"NICKE latdeg first {lat_deg.min()}, {lat_deg.max()}")

        lat_valid = (lat_deg <= 360) & (lat_deg >= -180)
        lon_valid = (lon_deg <= 360) & (lon_deg >= -180)

        lat_deg[lat_valid] = np.deg2rad(lat_deg[lat_valid])
        lon_deg[lon_valid] = np.deg2rad(lon_deg[lon_valid])

        print(f"NICKE latdeg second {lat_deg.min()}, {lat_deg.max()}")


        container.replace('gridLatitude', lat_deg, cat)
        container.replace('gridLongitude', lon_deg, cat)


    def _derive_stationIdentification(self, container, cat):
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
   
        container.add('stationIdentification', stid, said_paths, cat)
        #print("NICKE NEW CONTAINERLIST ", container.list()) 


    #def Co_imph(impp, elrc, geodu):
    def _derive_imph_and_bnda(self, container, cat):
        self.log.info("Get required variables")
        impp1 = container.get('impactParameterRO_roseq2repl1', cat).astype(np.float32)
        impp2 = container.get('impactParameterRO_roseq2repl2', cat).astype(np.float32)
        impp3 = container.get('impactParameterRO_roseq2repl3', cat).astype(np.float32)
        elrc = container.get('earthRadiusCurvature', cat)
        geodu = container.get('geoidUndulation', cat)

        mefr1 = container.get('frequency__roseq2repl1', cat)
        mefr2 = container.get('frequency__roseq2repl2', cat)
        mefr3 = container.get('frequency__roseq2repl3', cat)
        bnda1 = container.get('bendingAngle_roseq2repl1', cat)
        bnda2 = container.get('bendingAngle_roseq2repl2', cat)
        bnda3 = container.get('bendingAngle_roseq2repl3', cat)
        bndaoe1 = container.get('obsErrorBendingAngle1', cat)
        bndaoe2 = container.get('obsErrorBendingAngle2', cat)
        bndaoe3 = container.get('obsErrorBendingAngle3', cat)

        self.log.info("Calculate imph")
        imph1 = (impp1 - elrc - geodu).astype(np.float32)
        imph2 = (impp2 - elrc - geodu).astype(np.float32)
        imph3 = (impp3 - elrc - geodu).astype(np.float32)
   
        self.log.info("Overwrite values for mefr, bnda, impp, imph, bndaoe"
        for i in range(len(impp1)):
            if (mefr2[i] == 0.0):
                bnda1[i] = bnda2[i]
                mefr1[i] = mefr2[i]
                impp1[i] = impp2[i]
                imph1[i] = imph2[i]
                bndaoe1[i] = bndaoe2[i]
            if (mefr3[i] == 0.0):
                bnda1[i] = bnda3[i]
                mefr1[i] = mefr3[i]
                impp1[i] = impp3[i]
                imph1[i] = imph3[i]
                bndaoe1[i] = bndaoe3[i]       
    


        #container replace
        self.log.info("imph container replace")


        #container add
        self.log.info("imph container add")


        #container remove
        self.log.info("imph container remove")
        #container.remove('impactParameterRO_roseq2repl1')
        #container.remove('impactParameterRO_roseq2repl2')
        #container.remove('impactParameterRO_roseq2repl3')


        #return imph

        #satId = container.get('satelliteId', cat)
        #if not satId.size:
        #    self.log.warning(f'category {cat[0]} does not exist in input file')
        #    add_dummy_variable(container, 'pressure', cat, 'latitude')
        #    return
#
#        latitude = container.get('latitude', cat)
#        paths = container.get_paths('latitude', cat)
#
#        pressure = np.full_like(latitude, 0)
        #container.add('pressure', pressure, paths, cat)   

        

    #for i in range(len(degrees)):
    #    if degrees[i] <= 360 and degrees[i] >= -180:
    #        degrees[i] = np.deg2rad(degrees[i])
    #rad = degrees











#    def _add_wind_descriptions(self, description):
#        description.add_variables([
#            {
#                'name': 'ObsType/windEastward',
#                'source': 'obstype_uwind',
#                'units': '1',
#                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
#            },
#            {
#                'name': 'ObsType/windNorthward',
#                'source': 'obstype_vwind',
#                'units': '1',
#                'longName': 'Observation Type based on Satellite-derived Wind Computation Method and Spectral Band',
#            },
#            {
#                'name': 'ObsValue/windEastward',
#                'source': 'windEastward',
#                'units': 'm s-1',
#                'longName': 'Eastward Wind Component',
#            },
#            {
#                'name': 'ObsValue/windNorthward',
#                'source': 'windNorthward',
#                'units': 'm s-1',
#                'longName': 'Northward Wind Component',
#            }])
#
#    def _add_quality_info_and_gen_app_descriptions(self, description):
#        description.add_variables([
#            {
#                'name': 'MetaData/windGeneratingApplication',
#                'source': 'windGeneratingApplication',
#                'units': '1',
#                'longName': 'Wind Generating Application',
#            },
#            {
#                'name': 'MetaData/qualityInformationWithoutForecast',
#                'source': 'qualityInformationWithoutForecast',
#                'units': '1',
#                'longName': 'Quality Information Without Forecast',
#            }])
#
#
#    # Methods that are used to extend the obs data container
#    def _add_wind_obs(self, container, cat):
#        # Add new variables: ObsType/windEastward & ObsType/windNorthward
#        swcm = container.get('windComputationMethod', cat)
#        chanfreq = container.get('sensorCentralFrequency', cat)
#
#        if swcm.size == 0:
#            self.log.warning(f'category {cat[0]} does not exist in input file')
#            paths = container.get_paths('variables/windComputationMethod', cat)
#            obstype = container.get('variables/windComputationMethod', cat)
#            container.add('variables/obstype_uwind', obstype, paths, cat)
#            container.add('variables/obstype_vwind', obstype, paths, cat)
#
#            paths = container.get_paths('variables/windSpeed', cat)
#            wob = container.get('variables/windSpeed', cat)
#            container.add('variables/windEastward', wob, paths, cat)
#            container.add('variables/windNorthward', wob, paths, cat)
#            return
#
#        # self.log.debug(f'swcm min/max = {swcm.min()} {swcm.max()}')
#        self.log.debug('chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')
#
#        obstype = self._get_obs_type(swcm, chanfreq)
#
#        self.log.debug(f'obstype = {obstype}')
#        self.log.debug(f'obstype min/max =  {obstype.min()} {obstype.max()}')
#
#        paths = container.get_paths('windComputationMethod', cat)
#        container.add('obstype_uwind', obstype, paths, cat)
#        container.add('obstype_vwind', obstype, paths, cat)
#
#        # Add new variables: ObsValue/windEastward & ObsValue/windNorthward
#        wdir = container.get('windDirection', cat)
#        wspd = container.get('windSpeed', cat)
#
#        self.log.debug(f'wdir min/max = {wdir.min()} {wdir.max()}')
#        self.log.debug(f'wspd min/max = {wspd.min()} {wspd.max()}')
#
#        uob, vob = self._compute_wind_components(wdir, wspd)
#
#        self.log.debug(f'uob min/max = {uob.min()} {uob.max()}')
#        self.log.debug(f'vob min/max = {vob.min()} {vob.max()}')
#
#        paths = container.get_paths('windSpeed', cat)
#        container.add('windEastward', uob, paths, cat)
#        container.add('windNorthward', vob, paths, cat)
#
#    def _add_quality_info_and_gen_app(self, findQi, container, cat):
#        # Add new variables: MetaData/windGeneratingApplication and qualityInformationWithoutForecast
#        gnap2D = container.get('generatingApplication', cat)
#        pccf2D = container.get('qualityInformation', cat)
#        satId = container.get('satelliteId', cat)
#
#        if not satId.size:
#            paths = container.get_paths('windComputationMethod', cat)
#            dummy = container.get('windSpeed', cat)
#            container.add('windGeneratingApplication', dummy, paths, cat)
#            container.add('qualityInformationWithoutForecast', dummy, paths, cat)
#            return
#
#        gnap, qifn = self._get_quality_info_and_gen_app(findQi, gnap2D, pccf2D, satId)
#
#        self.log.debug(f'gnap min/max = {gnap.min()} {gnap.max()}')
#        self.log.debug(f'qifn min/max = {qifn.min()} {qifn.max()}')
#
#        paths = container.get_paths('windComputationMethod', cat)
#        container.add('windGeneratingApplication', gnap, paths, cat)
#        container.add('qualityInformationWithoutForecast', qifn, paths, cat)
#
#    def _get_obs_type(self, swcm, chan_freq=0):
#        """
#        Determine the observation type based on `swcm` and `chanfreq`.
#
#        Parameters:
#            swcm (array-like): Switch mode values.
#            chanfreq (array-like): Channel frequency values (Hz).
#
#        Returns:
#            numpy.ndarray: Observation type array.
#
#        Raises:
#            ValueError: If any `obstype` is unassigned.
#        """
#
#        raise NotImplementedError('Method _get_obs_type must be implemented in derived classes')
#
#    # Private methods
#    def _compute_wind_components(self, wdir, wspd):
#        """
#        Compute the U and V wind components from wind direction and wind speed.
#
#        Parameters:
#            wdir (array-like): Wind direction in degrees (meteorological convention: 0° = North, 90° = East).
#            wspd (array-like): Wind speed.
#
#        Returns:
#            tuple: U and V wind components as numpy arrays with dtype float32.
#        """
#        wdir_rad = np.radians(wdir)  # Convert degrees to radians
#        u = -wspd * np.sin(wdir_rad)
#        v = -wspd * np.cos(wdir_rad)
#
#        return u.astype(np.float32), v.astype(np.float32)
#
#    def _get_quality_info_and_gen_app(self, findQi, gnap2D, pccf2D):
#        # For NOAA VIIRS data, qi w/o forecast (qifn) is packaged in same
#        # vector of qi with ga = 5 (EUMETSAT QI without forecast). Must
#        # conduct a search and extract the correct vector for gnap and qi
#
#        # 1. Find dimension-sizes of ga and qi (should be the same!)
#        gDim1, gDim2 = np.shape(gnap2D)
#        qDim1, qDim2 = np.shape(pccf2D)
#        self.log.info('Generating Application and Quality Information SEARCH')
#        self.log.debug( f'Dimension size of GNAP ({gDim1},{gDim2})')
#        self.log.debug( f'Dimension size of PCCF ({qDim1},{qDim2})')
#
#        # 2. Initialize gnap and qifn as None, and search for dimension of
#        #    ga with values of 5. If the same column exists for qi, assign
#        #    gnap to ga[:,i] and qifn to qi[:,i], else raise warning that no
#        #    appropriate GNAP/PCCF combination was found
#        gnap = None
#        qifn = None
#        for i in range(gDim2):
#            if np.unique(gnap2D[:, i].squeeze()) == find_qi:
#                if i <= qDim2:
#                    self.log.info(f'GNAP/PCCF found for column {i}')
#                    gnap = gnap2D[:, i].squeeze()
#                    qifn = pccf2D[:, i].squeeze()
#                else:
#                    self.log.info(f'ERROR: GNAP column {i} outside of PCCF dimension {qDim2}')
#        if (gnap is None) & (qifn is None):
#            raise ValueError(f'GNAP == {findQI} NOT FOUND OR OUT OF PCCF DIMENSION-RANGE, WILL FAIL!')
#        # If EE is needed, key search on np.unique(gnap2D[:,i].squeeze()) == 7 instead
#        # NOTE: Make sure to return np.float32 or np.int32 types as appropriate!!!
#        return gnap.astype(np.int32), qifn.astype(np.int32)
