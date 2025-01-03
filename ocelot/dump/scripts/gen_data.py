import os
import sys
import argparse
from datetime import datetime, timedelta

import bufr
from bufr.encoders import zarr, netcdf
import emctank


OUTPUT_PATH='/scratch1/NCEPDEV/da/Ronald.McLaren/shared/ocelot/data'

def create_data_for_day(date:datetime, out_type:str):
    bufr.mpi.App(sys.argv)
    comm = bufr.mpi.Comm("world")

    start_datetime = date
    end_datetime = date + timedelta(hours=23, minutes=59, seconds=59)
    
    parameters = emctank.Parameters()
    parameters.start_time = start_datetime
    parameters.stop_time = end_datetime
    #parameters.category = ['n20']  # (Optional) Only keeps data from the specified category.
    parameters.gather = True   # (Optional) Gather data from all ranks to rank 0. 
                               # Framework preserves the order of the data.
    
    container, description = emctank.run(comm, "atms", parameters)

    if container is None:
        raise ValueError("No data found")

    if comm.rank() == 0:
        #file_exp = 'zarr' if out_type=='zarr' else 'nc'
        #output_path = os.path.join(OUTPUT_PATH, out_type, "atms", f'atms_{{splits/satId}}_{date.strftime("%Y_%m_%d")}.{file_exp}')
        
        #if out_type == "netcdf":
        #    netcdf.Encoder(description).encode(container, output_path)
        #elif out_type == "zarr":
        #    zarr.Encoder(description).encode(container, output_path)


        output_path = os.path.join(OUTPUT_PATH, "netcdf", "atms", f'atms_{{splits/satId}}_{date.strftime("%Y_%m_%d")}.nc')
        netcdf.Encoder(description).encode(container, output_path)
        print(f"Output written to {output_path}")
        
        output_path = os.path.join(OUTPUT_PATH, "zarr", "atms", f'atms_{{splits/satId}}_{date.strftime("%Y_%m_%d")}.zarr')
        zarr.Encoder(description).encode(container, output_path)
        print(f"Output written to {output_path}")

        sys.stdout.flush()

def create_data(start_date:datetime, end_date:datetime, out_type:str):
    date = start_date
    day = timedelta(days=1)

    while date <= end_date:
        create_data_for_day(date, out_type)
        date += day


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('start_date')
    parser.add_argument('end_date')
    parser.add_argument('-t', '--type',  choices=['zarr', 'netcdf'], default='zarr')

    args = parser.parse_args()

    start_date = datetime.strptime(args.start_date, "%Y-%m-%d")
    end_date = datetime.strptime(args.end_date, "%Y-%m-%d")

    create_data(start_date, end_date, args.type)

