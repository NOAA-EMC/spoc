import bufr
from pyioda.ioda.Engines.Bufr import Encoder
import numpy as np
import warnings
from wxflow import Logger

def get_description(mapping_path, update=False):

    description = bufr.encoders.Description(mapping_path)

    if update:

        description.add_variable(name='ObsType/windEastward',
                                 source='variables/obstype_uwind',
                                  units='1',
                                  longName='Observation Type based on Satellite-derived Wind Computation Method and Spectral Band')

        description.add_variable(name='ObsType/windNorthward',
                                 source='variables/obstype_vwind',
                                 units='1',
                                 longName='Observation Type based on Satellite-derived Wind Computation Method and Spectral Band')

        description.add_variable(name='ObsValue/windEastward',
                                 source='variables/windEastward',
                                 units='m/s',
                                 longName='Eastward Wind Component')

        description.add_variable(name='ObsValue/windNorthward',
                                 source='variables/windNorthward',
                                 units='m/s',
                                 longName='Northward Wind Component')

    return description

def get_all_keys(nested_dict):

    keys = []
    for key, value in nested_dict.items():
        keys.append(key)
        if isinstance(value, dict):
            keys.extend(get_all_keys(value))
    return keys

def Compute_WindComponents_from_WindDirection_and_WindSpeed(wdir, wspd):

    uob = (-wspd * np.sin(np.radians(wdir))).astype(np.float32)
    vob = (-wspd * np.cos(np.radians(wdir))).astype(np.float32)

    return uob, vob

def Get_ObsType(swcm, chanfreq):

    obstype = swcm.copy()

    # Use numpy vectorized operations
    obstype = np.where(swcm == 5, 247, obstype)  # WVCA/DL
    obstype = np.where(swcm == 3, 246, obstype)  # WVCT
    obstype = np.where(swcm == 2, 251, obstype)  # VIS
    obstype = np.where(swcm == 1, 245, obstype)  # IRLW

    condition = np.logical_and(swcm == 1, chanfreq >= 50000000000000.0)  # IRSW
    obstype = np.where(condition, 240, obstype)

    if not np.any(np.isin(obstype, [247, 246, 251, 245, 240])):
        raise ValueError("Error: Unassigned ObsType found ... ")

    return obstype

def create_obs_group(input_path, mapping_path, category, env, log_level=None ):

    # Set MPI parameters
    comm = bufr.mpi.Comm(env["comm_name"])
    rank = comm.rank() 

    # Set logger level
    if log_level is None:
        log_level = 'INFO'
    elif log_level not in ['DEBUG', 'INFO']:
        raise ValueError("Invalid log level! log_level must be 'DEBUG' or 'INFO'.")

    # Initialize logger
    logger = Logger('BUFR2IODA_satwnd_amv_goes.py', level=log_level, colored_log=False)

    # Generate keys for cache
    inputKey = input_path
    mappingKey = mapping_path

    if rank == 0: logger.info(f'Create obs group for category: {category}')

    # Check the cache for the data and return it if it exists
    if rank == 0: logger.debug(f'Check if bufr.DataCache exists? {bufr.DataCache.has(inputKey, mappingKey)}')
    if bufr.DataCache.has(inputKey, mappingKey):
        container = bufr.DataCache.get(inputKey, mappingKey)
        if rank == 0: logger.info(f'Encode {category} from cache')
        data = Encoder(get_description(mapping_path, update=True)).encode(container)[(category,)]
        if rank == 0: logger.info(f'Mark {category} as finished in the cache')
        bufr.DataCache.mark_finished(inputKey, mappingKey, [category])
        if rank == 0: logger.info(f'Return the encoded data for {category}')
        return data

    if rank == 0: logger.info(f'Get all cagegories info cache')
    # If cache does not exist, get data into cache 
    # Get container from mapping file first
    if rank == 0: logger.info(f'Get container from bufr parser')
    container = bufr.Parser(input_path, mapping_path).parse(comm)
    if rank == 0: logger.debug('container list (original):\n' + '\n'.join(str(item) for item in container.list()))
    if rank == 0: logger.debug(f'all_sub_categories =  {container.all_sub_categories()}')
    if rank == 0: logger.debug(f'category map =  {container.get_category_map()}')

    categories = container.all_sub_categories()
    for cat in categories:  

        if rank == 0: logger.debug(f'category = {cat}')
        description = get_description(mapping_path, update=True)

        satid = container.get('variables/satelliteId', cat)
        if satid.size == 0:
            if rank == 0: logger.info(f'category {cat[0]} does not exist in input file')
            paths = container.get_paths('variables/windComputationMethod', cat)
            obstype = container.get('variables/windComputationMethod', cat)
            container.add('variables/obstype_uwind', obstype, paths, cat)
            container.add('variables/obstype_vwind', obstype, paths, cat)

            paths = container.get_paths('variables/windSpeed', cat)
            wob = container.get('variables/windSpeed', cat)
            container.add('variables/windEastward', wob, paths, cat)
            container.add('variables/windNorthward', wob, paths, cat)

        else:
            # Add new variables: ObsType/windEastward & ObsType/windNorthward 
            swcm = container.get('variables/windComputationMethod', cat)
            chanfreq = container.get('variables/sensorCentralFrequency', cat)

            if rank == 0: logger.debug(f'swcm min/max = {swcm.min()} {swcm.max()}')
            if rank == 0: logger.debug(f'chanfreq min/max = {chanfreq.min()} {chanfreq.max()}')

            obstype = Get_ObsType(swcm, chanfreq)

            if rank == 0: logger.debug(f'obstype = {obstype}')
            if rank == 0: logger.debug(f'obstype min/max =  {obstype.min()} {obstype.max()}')

            paths = container.get_paths('variables/windComputationMethod', cat)
            container.add('variables/obstype_uwind', obstype, paths, cat)
            container.add('variables/obstype_vwind', obstype, paths, cat)

            # Add new variables: ObsValue/windEastward & ObsValue/windNorthward 
            wdir = container.get('variables/windDirection', cat)
            wspd = container.get('variables/windSpeed', cat)

            if rank == 0: logger.debug(f'wdir min/max = {wdir.min()} {wdir.max()}')
            if rank == 0: logger.debug(f'wspd min/max = {wspd.min()} {wspd.max()}')

            uob, vob = Compute_WindComponents_from_WindDirection_and_WindSpeed(wdir, wspd)

            if rank == 0: logger.debug(f'uob min/max = {uob.min()} {uob.max()}')
            if rank == 0: logger.debug(f'vob min/max = {vob.min()} {vob.max()}')

            paths = container.get_paths('variables/windSpeed', cat)
            container.add('variables/windEastward', uob, paths, cat)
            container.add('variables/windNorthward', vob, paths, cat)

    # Check
    if rank == 0: logger.debug('container list (updated):\n' + '\n'.join(str(item) for item in container.list()))
    if rank == 0: logger.debug(f'all_sub_categories {container.all_sub_categories()}')
     
    # Gather data from all tasks into all tasks. Each task will have the complete record 
    if rank == 0: logger.info(f'Gather data from all tasks into all tasks')
    container.all_gather(comm)

    if rank == 0: logger.info(f'Add container to cache')
    # Add the container to the cache
    bufr.DataCache.add(inputKey, mappingKey, container.all_sub_categories(), container)

    # Encode the data
    if rank == 0: logger.info(f'Encode {category}')
    data = Encoder(description).encode(container)[(category,)]

    if rank == 0: logger.info(f'Mark {category} as finished in the cache')
    # Mark the data as finished in the cache
    bufr.DataCache.mark_finished(inputKey, mappingKey, [category])

    if rank == 0: logger.info(f'Return the encoded data for {category}')
    return data
