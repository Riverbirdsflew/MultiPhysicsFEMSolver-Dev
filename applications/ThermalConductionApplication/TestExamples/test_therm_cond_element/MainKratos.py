import sys
import time
import importlib
from KratosMultiphysics.ThermalConductionApplication.thermal_conduction_analysis import ThermalConductionAnalysis
import KratosMultiphysics

if __name__ == "__main__":
    with open("basic_diffusion_test_stationary_parameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    global_model = KratosMultiphysics.Model()
    simulation = ThermalConductionAnalysis(global_model, parameters)
    simulation.Run()
