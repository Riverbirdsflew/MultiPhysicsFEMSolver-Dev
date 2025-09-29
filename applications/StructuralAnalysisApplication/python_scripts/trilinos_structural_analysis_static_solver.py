# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
from KratosMultiphysics.StructuralAnalysisApplication.trilinos_structural_analysis_solver import TrilinosStructuralSolver

def CreateSolver(model, custom_settings):
    return TrilinosStaticStructuralSolver(model, custom_settings)

class TrilinosStaticStructuralSolver(TrilinosStructuralSolver):
    """The trilinos structural mechanics static solver.

    For more information see:
    structural_analysis_solver.py
    trilinos_structural_analysis_solver.py
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosStaticStructuralSolver]:: ", "Construction finished")

    def _CreateScheme(self):
        return TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme()
