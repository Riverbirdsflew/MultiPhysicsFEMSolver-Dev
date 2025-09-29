import sys

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ThermalConductionApplication
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
from KratosMultiphysics.ThermalConductionApplication import python_solvers_wrapper_thermal_conduction           # Import the thermal Python solvers wrapper
from KratosMultiphysics.StructuralAnalysisApplication import python_solvers_wrapper_structural       # Import the structure Python solvers wrapper
from KratosMultiphysics.ConcentrationDiffusionApplication import python_solvers_wrapper_concentration_diffusion       # Import the Diffusion Python solvers wrapper

import KratosMultiphysics.ThermDiffStructPhaseApplication as KratosTDSP

def CreateSolver(main_model_part, custom_settings):
    return ThermDiffStructPhaseCoupledSolver(main_model_part, custom_settings)

class ThermDiffStructPhaseCoupledSolver(PythonSolver):
    
    @classmethod
    def GetDefaultParameters(cls):

        # Note that only the coupling settings are validated
        # The subdomain solver settings will be validated while instantiating these
        default_setting = KratosMultiphysics.Parameters("""
        {

            "parallel_type": "OpenMP",
            "start_time": 0.0,
            "end_time": 1.99e-1
            "solver_type" : "ThermDiffStructPhaseCoupled",
            "domain_size" : -1,
            "echo_level": 0,
            "calculation_flag_settings":{
                "phase_transformation_latent_calculation" : false,
                "thermal_strain_calculation": false,
                "phase_transformation_strain_calculation" : false,
                "trip_calculation" : false,
                "creep_calculation" : false
            },
            "structure_solver_settings": {
                "solver_type": "Static",
                "model_part_name"                 : "Structure",
                "domain_size"                     : 2,
                "echo_level"                      : 1,
                "analysis_type"                   : "non_linear",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                },
                "material_import_settings"        : {
                    "materials_filename" : "StructuralMaterials.json"
            }
            },
            "thermo_solver_settings":{
                "solver_type": "Transient",
                "analysis_type": "linear",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                }
            },
            "diffusion_solver_settings":{
                "solver_type": "Transient",
                "analysis_type": "linear",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "DiffuseMaterials.json"
                }
            },
            "time_integration_method": "implicit"
        }""")

        default_setting.AddMissingParameters(super().GetDefaultParameters())
        return default_setting
    
    def __init__(self, model, custom_settings):

        # Call the base Python solver constructor
        # Note that default settings in GetDefaultParameters() are validated in here
        super(ThermDiffStructPhaseCoupledSolver, self).__init__(model, custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()
        self.parallel_type = self.settings["parallel_type"].GetString()

        self.structural_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"],self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("ThermDiffStructPhaseCoupledSolver", "Structure solver construction finished.")

        self.thermal_solver = python_solvers_wrapper_thermal_conduction.CreateSolverByParameters(self.model,self.settings["thermo_solver_settings"],self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("ThermDiffStructPhaseCoupledSolver", "Thermo solver construction finished.")

        self.diffuse_solver = python_solvers_wrapper_concentration_diffusion.CreateSolverByParameters(self.model,self.settings["diffusion_solver_settings"],self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("ThermDiffStructPhaseCoupledSolver", "Diffusion solver construction finished.")

        KratosMultiphysics.Logger.PrintInfo("ThermDiffStructPhaseCoupledSolver", "Therm-Diff-Struct-phase solver construction finished.")


    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.structural_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        buffer_size_diffuse = self.diffuse_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal, buffer_size_diffuse)

    def AddVariables(self):
        ## Structure variables addition
        self.structure_solver.AddVariables()
        ## Thermp variables addition
        self.thermal_solver.AddVariables()
        # Diffusion solver variables addition
        self.diffuse_solver.AddVariables()

        KratosMultiphysics.MergeVariableListsUtility().Merge(self.structural_solver.main_model_part, self.thermal_solver.main_model_part)
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.structural_solver.main_model_part, self.diffuse_solver.main_model_part)
        
    def ImportModelPart(self):
        # structure solvers ImportModelPart() call
        self.structure_solver.ImportModelPart()

    def PrepareModelPart(self):
        self.structural_solver.PrepareModelPart()

        # Here the structural model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPartMixElemTherm(self.structural_solver.main_model_part,
                                        self.thermal_solver.main_model_part,
                                        "EulerianThermCond2D",
                                        "ThermalFace2D")
            modeler.GenerateModelPartMixElemDiff(self.structural_solver.main_model_part,
                                        self.diffuse_solver.main_model_part,
                                        "EulerianConcenDiff2D",
                                        "ConcentrationFace2D")
        elif self.domain_size == 3:
            modeler.GenerateModelPartMixElemTherm(self.structural_solver.main_model_part,
                                        self.thermal_solver.main_model_part,
                                        "EulerianThermCond3D",
                                        "ThermalFace3D")
            modeler.GenerateModelPartMixElemDiff(self.structural_solver.main_model_part,
                                        self.diffuse_solver.main_model_part,
                                        "EulerianConcenDiff3D",
                                        "ConcentrationFace3D")
        self.thermal_solver.PrepareModelPart()
        self.diffuse_solver.PrepareModelPart()

    def AddDofs(self):
        # Add DOFs
        self.thermal_solver.AddDofs()
        self.diffuse_solver.AddDofs()
        self.structural_solver.AddDofs()

    def Initialize(self):
        self.thermal_solver.Initialize()
        self.diffuse_solver.Initialize()
        self.structural_solver.Initialize()
        # Set Clalulation flags
        self._SetCalculationFlags(self.settings["calculation_flag_settings"])

    def AdvanceInTime(self, current_time):
        new_time = self.structural_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        # Initialize solution step of coupling solvers
        self.thermal_solver.InitializeSolutionStep()
        self.diffuse_solver.InitializeSolutionStep()
        self.structural_solver.InitializeSolutionStep()

    def Predict(self):
        pass

    def GetComputingModelPart(self):
        return self.structural_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def SaveRestart(self):
        pass

    def SolveSolutionStep(self):
        
        self.thermal_solver.InitializeSolutionStep()
        self.thermal_solver.Predict()
        KratosMultiphysics.Logger.PrintInfo("\t" + "Solving THERMAL part...")
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()
        
        self.diffuse_solver.InitializeSolutionStep()
        self.diffuse_solver.Predict()
        KratosMultiphysics.Logger.PrintInfo("\t" + "Solving DIFFUSE part...")
        diffuse_is_converged = self.diffuse_solver.SolveSolutionStep()

        self.structural_solver.InitializeSolutionStep()
        self.structural_solver.Predict()
        KratosMultiphysics.Logger.PrintInfo("\t" + "Solving STRUCTURAL part...")
        solid_is_converged = self.structural_solver.SolveSolutionStep()

        is_converged = solid_is_converged and thermal_is_converged and diffuse_is_converged

        return is_converged

    def FinalizeSolutionStep(self):
        self.thermal_solver.FinalizeSolutionStep()
        self.diffuse_solver.FinalizeSolutionStep()
        self.structural_solver.FinalizeSolutionStep()

    def Finalize(self):
        self.thermal_solver.Finalize()
        self.diffuse_solver.Finalize()
        self.structural_solver.Finalize()

    def SetEchoLevel(self, level):
        (self.thermal_solver).SetEchoLevel(level)
        (self.diffuse_solver).SetEchoLevel(level)
        (self.structural_solver).SetEchoLevel(level)

    def Clear(self):
        (self.structural_solver).Clear()
        (self.thermal_solver).Clear()
        (self.diffuse_solver).Clear()

    def Check(self):
        (self.structural_solver).Check()
        (self.thermal_solver).Check()
        (self.diffuse_solver).Check()

    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################
    
    def _SetCalculationFlags(self, custom_calculation_settings):
        """Set calculation flags for the solvers"""
        # Set the calculation flags for the solver
        calc_settins = custom_calculation_settings
        self.options = KratosMultiphysics.Flags()

        if( calc_settins["phase_transformation_strain_calculation"].GetBool() ):
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_PHASE_TRANSFORMATION_STRAIN, True)
        else:
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_PHASE_TRANSFORMATION_STRAIN, False)
        if( calc_settins["thermal_strain_calculation"].GetBool() ):
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_THERMAL_STRAIN, True)
        else:
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_THERMAL_STRAIN, False)
        if( calc_settins["trip_calculation"].GetBool() ):
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_TRIP, True)
        else:
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_TRIP, False)
        if( calc_settins["creep_calculation"].GetBool() ):
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CREEP, True)
        else:
            self.structural_solver.ProcessInfo.SetValue(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CREEP, False)
        if( calc_settins["phase_transformation_latent_calculation"].GetBool() ):
            self.thermal_solver.ProcessInfo.SetValue(KratosMultiphysics.ThermalConductionApplication.COMPUTE_PHASE_TRANSFORMATION_LATENT, True)
        else:
            self.thermal_solver.ProcessInfo.SetValue(KratosMultiphysics.ThermalConductionApplication.COMPUTE_PHASE_TRANSFORMATION_LATENT, False)