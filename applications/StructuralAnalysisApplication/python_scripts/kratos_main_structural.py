import KratosMultiphysics
from KratosMultiphysics.StructuralAnalysisApplication.structural_analysis import StructuralAnalysis

"""
For user-scripting it is intended that a new class is derived
from StructuralAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralAnalysis(model,parameters)
    simulation.Run()
