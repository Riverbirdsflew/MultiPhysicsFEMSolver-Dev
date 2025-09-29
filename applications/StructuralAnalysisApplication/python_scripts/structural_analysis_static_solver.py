# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralAnalysisApplication.structural_analysis_solver import StructuralSolver

def CreateSolver(model, custom_settings):
    return StaticStructuralSolver(model, custom_settings)

class StaticStructuralSolver(StructuralSolver):
    """The structural mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    See structural_analysis_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticStructuralSolver]:: ", "Construction finished")

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
