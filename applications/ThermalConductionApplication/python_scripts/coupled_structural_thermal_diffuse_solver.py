import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralAnalysisApplication as KratosSAA
import KratosMultiphysics.ThermalConductionApplication as ThermCond
import KratosMultiphysics.ConcentrationDiffusionApplication as ConcenDiff


# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledThermoDiffuseMechanicalSolver(main_model_part, custom_settings)

class CoupledThermoDiffuseMechanicalSolver(PythonSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            
                
                "solver_settings" : {
                "solver_type" : "ThermoDiffuseMechanicallyCoupled",
                "domain_size" : -1,
                "echo_level": 0,
                "structural_solver_settings": {
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
                "thermal_solver_settings": {
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "use_input_model_part"
                    },
                    "material_import_settings": {
                            "materials_filename": "ThermalMaterials.json"
                    }
                },
                "diffuse_solver_settings": {
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
            }
        }
        """)

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):

        super(CoupledThermoDiffuseMechanicalSolver, self).__init__(model, custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.StructuralAnalysisApplication import python_solvers_wrapper_structural
        self.structural_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structural_solver_settings"],"OpenMP")

        from KratosMultiphysics.ThermalConductionApplication import python_solvers_wrapper_thermal_conduction
        self.thermal_solver = python_solvers_wrapper_thermal_conduction.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")
        
        from KratosMultiphysics.ConcentrationDiffusionApplication import python_solvers_wrapper_concentration_diffusion
        self.diffuse_solver = python_solvers_wrapper_concentration_diffusion.CreateSolverByParameters(self.model,self.settings["diffuse_solver_settings"],"OpenMP")
        
        solver_type = self.settings["structural_solver_settings"]["solver_type"].GetString()
        self.is_dynamic = solver_type == "dynamic"

    def AddVariables(self):
        # Import the structural and thermal solver variables. Then merge them to have them in both structural and thermal solvers.
        self.structural_solver.AddVariables()
        self.thermal_solver.AddVariables()
        self.diffuse_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.structural_solver.main_model_part, self.thermal_solver.main_model_part)
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.structural_solver.main_model_part, self.diffuse_solver.main_model_part)
        #KratosMultiphysics.MergeVariableListsUtility().Merge(self.thermal_solver.main_model_part, self.diffuse_solver.main_model_part)

    def ImportModelPart(self):
        # Call the structural solver to import the model part from the mdpa
        self.structural_solver.ImportModelPart()

    def PrepareModelPart(self):
        self.structural_solver.PrepareModelPart()

        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        concentration_diffusion_settings = self.diffuse_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONCENTRATION_DIFFUSION_SETTINGS)

        # Here the structural model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "Element2D",
                                      "ThermalFace2D2N")
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      self.diffuse_solver.main_model_part,
                                      "Element2D",
                                      "ConcentrationFace2D2N")
            
            thermal_replace_settings = KratosMultiphysics.Parameters("""
                {
                    "element_name": "",
                    "condition_name": "ThermalFace2D2N"
                }
                """)
            diffuse_replace_settings = KratosMultiphysics.Parameters("""
                {
                    "element_name": "",
                    "condition_name": "ConcentrationFace2D2N"
                }
                """)
            
            for elem in self.thermal_solver.main_model_part.Elements:
                num_nodes = len(elem.GetNodes())
                if elem.Info == 'Element2D3N':  # todo: 用element.info()判断单元类型，用节点数可能产生冲突,对应单元名称加前缀或者后缀
                    thermal_replace_settings["element_name"].SetString("EulerianThermCond2D3N")
                elif elem.Info == 'Element2D4N':
                    thermal_replace_settings["element_name"].SetString("EulerianThermCond2D4N")
                
                # 执行替换
                KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.thermal_solver.main_model_part, thermal_replace_settings).Execute()
            
            for elem in self.diffuse_solver.main_model_part.Elements:
                num_nodes = len(elem.GetNodes())
                if elem.Info == 'Element2D3N':
                    diffuse_replace_settings["element_name"].SetString("EulerianConcenDiff2D3N")
                elif elem.Info == 'Element2D4N':
                    diffuse_replace_settings["element_name"].SetString("EulerianConcenDiff2D4N")
                # 执行替换
                KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.diffuse_solver.main_model_part, diffuse_replace_settings).Execute()
            #KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.thermal_solver.main_model_part,self.thermal_solver._get_element_condition_replace_settings()).Execute()
            #KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.diffuse_solver.main_model_part,self.diffuse_solver._get_element_condition_replace_settings()).Execute()
        else:
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianThermCondElement3D8N",
                                      "ThermalFace3D4N")
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      self.diffuse_solver.main_model_part,
                                      "EulerianConcenDiff3D8NElement",
                                      "ConcentrationFace3D4N")
            thermal_replace_settings = KratosMultiphysics.Parameters("""
                {
                    "element_name": "",
                    "condition_name": "ThermalFace3D3N"
                }
                """)
            diffuse_replace_settings = KratosMultiphysics.Parameters("""
                {
                    "element_name": "",
                    "condition_name": "ConcentrationFace3D3N"
                }
                """)

            # for elem in self.thermal_solver.main_model_part.Elements:
            #     num_nodes = len(elem.GetNodes())
            #     if elem.Info == 'Element3D4N':
            #         thermal_replace_settings["element_name"].SetString("EulerianThermCond3D4N")
            #     elif elem.Info == 'Element3D8N':
            #         thermal_replace_settings["element_name"].SetString("EulerianThermCond3D8N")
            #         thermal_replace_settings["condition_name"].SetString("ThermalFace3D4N")
            #     KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.thermal_solver.main_model_part, thermal_replace_settings).Execute()

            # for elem in self.diffuse_solver.main_model_part.Elements:
            #     num_nodes = len(elem.GetNodes())
            #     if elem.Info == 'Element3D4N':
            #         diffuse_replace_settings["element_name"].SetString("EulerianConcenDiff3D4N")
            #     elif elem.Info == 'Element3D8N':
            #         diffuse_replace_settings["element_name"].SetString("EulerianConcenDiff3D8N")
            #         diffuse_replace_settings["condition_name"].SetString("ConcentrationFace3D4N")
            #     KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.diffuse_solver.main_model_part, diffuse_replace_settings).Execute()

        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
        self.diffuse_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONCENTRATION_DIFFUSION_SETTINGS, concentration_diffusion_settings)

        self.thermal_solver.PrepareModelPart()
        self.diffuse_solver.PrepareModelPart()

    def AddDofs(self):
        self.structural_solver.AddDofs()
        self.thermal_solver.AddDofs()
        self.diffuse_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.structural_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.structural_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        buffer_size_diffuse = self.diffuse_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal, buffer_size_diffuse)

    def Initialize(self):
        self.structural_solver.Initialize()
        self.thermal_solver.Initialize()
        self.diffuse_solver.Initialize()

    def Clear(self):
        (self.structural_solver).Clear()
        (self.thermal_solver).Clear()
        (self.diffuse_solver).Clear()

    def Check(self):
        (self.structural_solver).Check()
        (self.thermal_solver).Check()
        (self.diffuse_solver).Check()

    def SetEchoLevel(self, level):
        (self.structural_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)
        (self.diffuse_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.structural_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
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

        self.RemoveConvectiveVelocity()

        return solid_is_converged and thermal_is_converged and diffuse_is_converged

    def FinalizeSolutionStep(self):
        self.structural_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()
        self.diffuse_solver.FinalizeSolutionStep()

    def RemoveConvectiveVelocity(self):
        if self.is_dynamic:
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VELOCITY, KratosMultiphysics.MESH_VELOCITY, self.thermal_solver.GetComputingModelPart(), self.thermal_solver.GetComputingModelPart(), 0)
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VELOCITY, KratosMultiphysics.MESH_VELOCITY, self.diffuse_solver.GetComputingModelPart(), self.diffuse_solver.GetComputingModelPart(), 0)
            
    def _get_thermal_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")

        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)

        ## Elements
        ## Note that we check for the elements that require substitution to allow for custom elements
        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        element_list = ["EulerianThermCond","LaplacianElement","MixedLaplacianElement","AdjointHeatDiffusionElement","QSConvectionDiffusionExplicit","DConvectionDiffusionExplicit","AxisymmetricEulerianThermalConduction"]
        if element_name in element_list:
            num_nodes_elements = 0
            if (len(self.main_model_part.Elements) > 0):
                for elem in self.main_model_part.Elements:
                    num_nodes_elements = len(elem.GetNodes())
                    break

            num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
            if not num_nodes_elements:
                num_nodes_elements = domain_size + 1

            name_string = f"{element_name}{domain_size}D{num_nodes_elements}N"
            self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        ## Conditions
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
        condition_list = ["FluxCondition","ThermalFace","AxisymmetricThermalFace","LineCondition","SurfaceCondition"]
        if condition_name in condition_list:
            num_nodes_conditions = 0
            if (len(self.main_model_part.Conditions) > 0):
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break

            num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
            if not num_nodes_conditions:
                num_nodes_conditions = domain_size

            name_string = f"{condition_name}{domain_size}D{num_nodes_conditions}N"
            self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]
        
    def _get_diffuse_element_condition_replace_settings(self):
        ## Get and check domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")

        ## Validate the replace settings
        default_replace_settings = self.GetDefaultParameters()["element_replace_settings"]
        self.settings["element_replace_settings"].ValidateAndAssignDefaults(default_replace_settings)

        ## Elements
        ## Note that we check for the elements that require substitution to allow for custom elements
        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        element_list = ["EulerianConcenDiff","LaplacianElement","MixedLaplacianElement","AdjointHeatDiffusionElement","QSConvectionDiffusionExplicit","DConvectionDiffusionExplicit","AxisymmetricEulerianConcentrationDiffusion"]
        if element_name in element_list:
            num_nodes_elements = 0
            if (len(self.main_model_part.Elements) > 0):
                for elem in self.main_model_part.Elements:
                    num_nodes_elements = len(elem.GetNodes())
                    break

            num_nodes_elements = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
            if not num_nodes_elements:
                num_nodes_elements = domain_size + 1

            name_string = f"{element_name}{domain_size}D{num_nodes_elements}N"
            self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        ## Conditions
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
        condition_list = ["ConcentrationFluxCondition","ConcentrationFace2D2N","AxisymmetricConcentrationFace","LineCondition","SurfaceCondition"]
        if condition_name in condition_list:
            num_nodes_conditions = 0
            if (len(self.main_model_part.Conditions) > 0):
                for cond in self.main_model_part.Conditions:
                    num_nodes_conditions = len(cond.GetNodes())
                    break

            num_nodes_conditions = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
            if not num_nodes_conditions:
                num_nodes_conditions = domain_size

            name_string = f"{condition_name}{domain_size}D{num_nodes_conditions}N"
            self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        return self.settings["element_replace_settings"]

