#!/usr/bin/env python

from pycube import CubeModel


class Model:
    """
    ATOM Cube Turbulence Model wrapper.
    Drives the C++ cCubeModel via the pycube Cython extension.
    """
    def __init__(self, config="cli/config_cube.xml"):
        self.cube = CubeModel()
        self.cube.load_config(config)
        self.config_xml = config
        self.time_slice = 0

    def print_config(self):
        print("\n\n\n\n   ATOM Cube Turbulence Model")
        print(f"\n   configuration file       : {self.config_xml}")
        print(f"   output path              : {self.cube.output_path}")

    def run(self, time_slice=0):
        self.time_slice = time_slice
        print(f"\n   Running turbulent flow around a cube  —  time-slice Ma = {self.time_slice}")
        self.cube.run()
        print(f"\n   Successfully terminated  —  time-slice Ma = {self.time_slice}\n")


if __name__ == "__main__":
    model = Model()
    model.print_config()
    model.run(0)
