
TEST_CASES = [
    {
        'map': 'bufr_atms',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.atms.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_atms_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_atms_npp.nc',
        'result': 'testrun/2022050100-bufr_atms_npp.nc'
    },
    {
        'map': 'bufr_cris-fsr',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.crisf4.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_cris_fsr_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_cris_fsr_npp.nc',
        'result': 'testrun/2022050100-bufr_cris_fsr_npp.nc'
    },
    {
        'map': 'bufr_iasi',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.mtiasi.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_iasi_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_iasi_metop-b.nc',
        'result': 'testrun/2022050100-bufr_iasi_metop-b.nc'
    },
    {
        'map': 'bufr_satwnd_amv_abi',
        'args' : {
            'input': 'testdata/2022050100-gdas.t00z.satwnd.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_satwnd_amv_abi_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_satwnd_amv_abi_goes-16.nc',
        'result': 'testrun/2022050100-bufr_satwnd_amv_abi_goes-16.nc'
    },
    {
        'map': 'bufr_satwnd_amv_ahi',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.satwnd.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_satwnd_amv_ahi_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_satwnd_amv_ahi_h9.nc',
        'result': 'testrun/2022050100-bufr_satwnd_amv_ahi_h9.nc'
    },
    {
        'map': 'bufr_satwnd_amv_avhrr',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.satwnd.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_satwnd_amv_avhrr_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_satwnd_amv_avhrr_n19.nc',
        'result': 'testrun/2022050100-bufr_satwnd_amv_avhrr_n19.nc'
    },
    {
        'map': 'bufr_satwnd_amv_leogeo',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.satwnd.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_satwnd_amv_leogeo_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_satwnd_amv_leogeo_multi.nc',
        'result': 'testrun/2022050100-bufr_satwnd_amv_leogeo_multi.nc'
    },
    {
        'map': 'bufr_satwnd_amv_modis',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.satwnd.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_satwnd_amv_modis_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_satwnd_amv_modis_terra.nc',
        'result': 'testrun/2022050100-bufr_satwnd_amv_modis_terra.nc'
    },
    {
        'map': 'bufr_satwnd_amv_seviri',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.satwnd.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_satwnd_amv_seviri_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_satwnd_amv_seviri_m10.nc',
        'result': 'testrun/2022050100-bufr_satwnd_amv_seviri_m10.nc'
    },
    {
        'map': 'bufr_scatwnd_ascat',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.ssmisu.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_scatwnd_ascat_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_scatwnd_ascat_metop-b.nc',
        'result': 'testrun/2022050100-bufr_scatwnd_ascat_metop-b.nc'
    },
    {
        'map': 'bufr_sfcsno',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.ssmisu.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_sfcsno.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_sfcsno.nc',
        'result': 'testrun/2022050100-bufr_sfcsno.nc'
    },
    {
        'map': 'bufr_ssmis',
        'args': {
            'input': 'testdata/2022050100-gdas.t00z.ssmisu.tm00.bufr_d',
            'output': 'testrun/2022050100-bufr_ssmis_{splits/satId}.nc'
        },
        'cmp': 'testoutput/2022050100-bufr_ssmis_f16.nc',
        'result': 'testrun/2022050100-bufr_ssmis_f16.nc'
    },
    {
        'map': 'prepbufr_adpsfc',
        'args': {
            'input': 'testdata/ 2022050100-gdas.t00z.prepbufr',
            'output': 'testrun/2022050100-prepbufr_adpsfc.nc'
        },
        'cmp': 'testoutput/2022050100-prepbufr_adpsfc.nc',
        'result': 'testrun/2022050100-prepbufr_adpsfc.nc'
    },
    {
        'map': 'prepbufr_sfcshp',
        'args': {
            'input': 'testdata/ 2022050100-gdas.t00z.prepbufr',
            'output': 'testrun/2022050100-prepbufr_sfcshp.nc'
        },
        'cmp': 'testoutput/2022050100-prepbufr_sfcshp.nc',
        'result': 'testrun/2022050100-prepbufr_sfcshp.nc'
    },
]
