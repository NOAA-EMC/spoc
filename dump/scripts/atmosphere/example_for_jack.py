#!/usr/bin/env python3
import bufr
from bufr.encoders import netcdf
import numpy as np

def main():
    container = bufr.DataContainer()
    sin = np.sin(np.linspace(0, 2 * np.pi, 100))
    container.add('sin', sin, ['*'])
    yamlpath = './example_for_jack.yaml'
    description = bufr.encoders.Description(yamlpath)
    netcdf.Encoder(description).encode(container, './out.nc')

if __name__ == "__main__":
    main()
