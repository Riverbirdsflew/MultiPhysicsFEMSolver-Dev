import KratosMultiphysics
# Import base class file
from KratosMultiphysics.ThermalConductionApplication import thermal_conduction_transient_solver
from KratosMultiphysics.ThermalConductionApplication import phase_transformation_model

def CreateSolver(main_model_part, custom_settings):
    return ThermalPhaseFieldSolver(main_model_part, custom_settings)


class ThermalPhaseFieldSolver(thermal_conduction_transient_solver.ThermalConductionTransientSolver):
    def Initialize(self):
        super().Initialize()
        self.phase_model = PhaseTransformationModel(self.model_part)
        print("Thermal-phase solver initialized.")
    
    def SolveSolutionStep(self):
        super().InitializeSolutionStep()
        super().Predict()

        KratosMultiphysics.Logger.PrintInfo("\t" + "Solving THERMAL part...")
        thermal_is_converged = super().SolveSolutionStep()
        self.phase_model.ComputePhaseFractions()
        
        return thermal_is_converged
    

