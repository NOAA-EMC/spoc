import json
import yaml
import os
import sys
from collections import OrderedDict


def bufr_filename(cycle_datetime, cycle_type, hh, data_format):
    return f"{cycle_datetime}-{cycle_type}.t{hh}z.{data_format}.tm00.bufr_d"


def ioda_filename(cycle_type, hh, descriptor, cycle_datetime):
    return f"{cycle_type}.t{hh}z.insitu_{descriptor}.{cycle_datetime}.nc4"



# configuration file can be either json or yaml
# this config class provides the functions that determine
# the names and the paths of the bufr input and the ioda output files
# these functions can be overridden in the converter


# Custom YAML dumper to preserve dictionary order
class OrderedDumper(yaml.SafeDumper):
    pass

def _dict_representer(dumper, data):
    return dumper.represent_mapping(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        data.items()
    )


class Bufr2iodaConfig:
    # def __init__(self, script_name, config_file, platform_description):
        # self.script_name = script_name
        # self.platform_description = platform_description
        # # read_config_file(config_file)

    def __init__(self):
        OrderedDumper.add_representer(OrderedDict, _dict_representer)

    def read_config_file(self, config_file):
        _, file_extension = os.path.splitext(config_file)
        if file_extension == ".json":
            with open(config_file, "r") as file:
                config = json.load(file)
            self.set(config)
        elif file_extension == ".yaml":
            with open(config_file, "r") as file:
                config = yaml.safe_load(file)
            self.set(config)
        else:
            print("Fatal error: Unknown file extension = ", file_extension)
            sys.exit(1)

    def set(self, config):
        # Get parameters from configuration
        self.data_format = config["data_format"]
        self.source = config["source"]
        self.data_type = config["data_type"]
        self.data_description = config["data_description"]
        self.data_provider = config["data_provider"]
        self.cycle_type = config["cycle_type"]
        self.cycle_datetime = config["cycle_datetime"]
        self.dump_dir = config["dump_directory"]
        self.ioda_dir = config["ioda_directory"]
        self.ocean_basin = config["ocean_basin"]

        self.yyyymmdd = self.cycle_datetime[0:8]
        self.hh = self.cycle_datetime[8:10]

        # General Information
        self.converter = 'BUFR to IODA Converter'

    def ocean_basin_nc_file_path(self):
        return self.ocean_basin

    def bufr_filename(self):
        return f"{self.cycle_datetime}-{self.cycle_type}.t{self.hh}z.{self.data_format}.tm00.bufr_d"

    def bufr_filepath(self):
        return os.path.join(self.dump_dir, self.bufr_filename())

    def ioda_filename(self, descriptor):
        return f"{self.cycle_type}.t{self.hh}z.insitu_{descriptor}.{self.cycle_datetime}.nc4"

    def ioda_filepath(self, descriptor):
        return os.path.join(self.ioda_dir, self.ioda_filename(descriptor))

    def create_config_file(self, data_format, subsets, 
                            data_type, data_description, 
                            cycle_type, cycle_datetime, 
                            dump_dir, ioda_dir, 
                            ocean_basin, output_path):
        """
        Create a YAML config file for BUFR to IODA conversion with keys in specified order.
        
        Args:
            data_format (str): Data format (e.g., 'bathy').
            subsets (str): BUFR subsets (e.g., 'BATHY').
            data_type (str): Data type (e.g., 'bathy').
            cycle_type (str): Cycle type (e.g., 'gdas').
            cycle_datetime (str): Cycle date-time (e.g., '2021063006').
            dump_dir (str): Path to BUFR input directory.
            ioda_dir (str): Path to IODA output directory.
            ocean_basin (str): Path to ocean basin file.
            output_path (str): Path for output YAML file.
        
        Raises:
            ValueError: If string inputs are empty or invalid.
            FileNotFoundError: If paths are invalid.
            OSError: If writing output file fails.
        """
        # Validate string inputs
        for param, value in [
            ('data_format', data_format), ('subsets', subsets),
            ('data_type', data_type), ('cycle_type', cycle_type),
            ('cycle_datetime', cycle_datetime),
            ('data_description', data_description)
        ]:
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{param} must be a non-empty string")

        # Validate paths
        if not os.path.isdir(dump_dir):
            raise FileNotFoundError(f"BUFR directory not found: {dump_dir}")
        if not os.path.isdir(ioda_dir):
            raise FileNotFoundError(f"IODA directory not found: {ioda_dir}")
        if not os.path.isfile(ocean_basin):
            raise FileNotFoundError(f"Ocean basin file not found: {ocean_basin}")

        # Create config with OrderedDict to preserve key order
        config = OrderedDict([
            ('data_format', data_format),
            ('subsets', subsets),
            ('source', 'NCEP data tank'),
            ('data_type', data_type),
            ('cycle_type', cycle_type),
            ('cycle_datetime', cycle_datetime),
            ('dump_directory', dump_dir),
            ('ioda_directory', ioda_dir),
            ('ocean_basin', ocean_basin),
            ('data_description', data_description),
            ('data_provider', 'U.S. NOAA')
        ])

        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Write YAML file with 2-space indentation, preserving order
        with open(output_path, 'w') as f:
            yaml.dump(config, f, Dumper=OrderedDumper, default_flow_style=False, indent=2)

        # print(f"Created config file: {output_path}")



    def replace_config_placeholders(config_path, bufr_dir, ioda_dir, ocean_basin_file, output_path):
        """
        Replace placeholders in a YAML config file and save to output_path.
        
        Args:
            config_path (str): Path to input YAML config file.
            bufr_dir (str): Path to BUFR input directory (replaces __BUFRINPUTDIR__).
            ioda_dir (str): Path to IODA output directory (replaces __IODAOUTPUTDIR__).
            ocean_basin_file (str): Path to ocean basin file (replaces __OCEANBASIN__).
            output_path (str): Path for output YAML file.
        
        Raises:
            FileNotFoundError: If config_path or directories are invalid.
            yaml.YAMLError: If YAML parsing fails.
            OSError: If writing output file fails.
        """
        # Validate input paths
        if not os.path.isfile(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")
        if not os.path.isdir(bufr_dir):
            raise FileNotFoundError(f"BUFR directory not found: {bufr_dir}")
        if not os.path.isdir(ioda_dir):
            raise FileNotFoundError(f"IODA directory not found: {ioda_dir}")
        if not os.path.isfile(ocean_basin_file):
            raise FileNotFoundError(f"Ocean basin file not found: {ocean_basin_file}")

        # Read YAML config
        with open(config_path, 'r') as f:
            config_text = f.read()

        # Replace placeholders
        config_text = config_text.replace('__BUFRINPUTDIR__', bufr_dir)
        config_text = config_text.replace('__IODAOUTPUTDIR__', ioda_dir)
        config_text = config_text.replace('__OCEANBASIN__', ocean_basin_file)

        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Write new YAML file
        with open(output_path, 'w') as f:
            f.write(config_text)

        # print(f"Created config file: {output_path}")



if __name__ == "__main__":

    # test the above function:
    config_path = "/work/noaa/da/edwardg/spoc/test/testconfig/bufr2ioda_insitu_profile_bathy_2021063006.yaml.in"
    bufr_dir = "/work/noaa/da/edwardg/"
    ioda_dir = "/work/noaa/da/edwardg/"
    ocean_basin_file = "/work/noaa/da/edwardg/spoc/test/jocean_basin.txt"
    output_path = "/work/noaa/da/edwardg/spoc/test/testconfig/jbufr2ioda_insitu_profile_bathy_2021063006.yaml"

    replace_config_placeholders(config_path, bufr_dir, ioda_dir, ocean_basin_file, output_path)



    # test the other function
    try:
        create_config_file(
            data_format="bathy",
            subsets="BATHY",
            data_type="bathy",
            cycle_type="gdas",
            cycle_datetime="2021063006",
            dump_dir="/work/noaa/da/edwardg",
            ioda_dir="/work/noaa/da/edwardg",
            ocean_basin="/work/noaa/da/edwardg/spoc/test/jocean_basin.txt",
            output_path="/work/noaa/da/edwardg/spoc/test/testconfig/jjbufr2ioda_insitu_profile_bathy_2021063006.yaml"
        )
    except Exception as e:
        print(f"Error: {e}")
