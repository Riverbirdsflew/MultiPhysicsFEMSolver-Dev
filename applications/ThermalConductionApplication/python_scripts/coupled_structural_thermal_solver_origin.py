import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralAnalysisApplication as KratosSAA
import KratosMultiphysics.ThermalConductionApplication as ThermCond

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledThermoMechanicalSolver(main_model_part, custom_settings)

class CoupledThermoMechanicalSolver(PythonSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermoMechanicallyCoupled",
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
            "time_integration_method": "implicit"
        }
        """)

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):

        super(CoupledThermoMechanicalSolver, self).__init__(model, custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.StructuralAnalysisApplication import python_solvers_wrapper_structural
        self.structural_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structural_solver_settings"],"OpenMP")

        from KratosMultiphysics.ThermalConductionApplication import python_solvers_wrapper_thermal_conduction
        self.thermal_solver = python_solvers_wrapper_thermal_conduction.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")
        solver_type = self.settings["structural_solver_settings"]["solver_type"].GetString()
        self.is_dynamic = solver_type == "dynamic"

    def AddVariables(self):
        # Import the structural and thermal solver variables. Then merge them to have them in both structural and thermal solvers.
        self.structural_solver.AddVariables()
        self.thermal_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.structural_solver.main_model_part, self.thermal_solver.main_model_part)

    def ImportModelPart(self):
        # Call the structural solver to import the model part from the mdpa
        self.structural_solver.ImportModelPart()

    def PrepareModelPart(self):
        self.structural_solver.PrepareModelPart()
        
        #for elem in self.structural_solver.main_model_part.Elements:
            #print(f"Element ID: {elem.Id}, Number of Nodes: {elem}")

        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        
        # 定义替换设置参数
        replace_settings = KratosMultiphysics.Parameters("""
        {
            "element_name": "",
            "condition_name": "ThermalFace2D2N"
        }
        """)

        # Here the structural model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                        self.thermal_solver.main_model_part,
                                        "EulerianThermCondElement"
                                        "ThermalFace2D2N")

            # 通过单元节点数来分别替换三节点和四节点单元
            for elem in self.thermal_solver.main_model_part.Elements:
                num_nodes = len(elem.GetNodes())
                if elem.Info() == "Element2D3N":
                    replace_settings["element_name"].SetString("EulerianThermCondElement2D3N")
                elif elem.Info() == "Element2D4N":
                    replace_settings["element_name"].SetString("EulerianThermCondElement2D4N")
                
                # 执行替换
                KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.thermal_solver.main_model_part, replace_settings).Execute()
                print(f"Element ID: {elem.Id}, Number of Nodes: {elem}")
            
        #else:
            #modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      #self.thermal_solver.main_model_part,
                                      #"EulerianThermCond3D",
                                      #"ThermalFace3D3N")
            #KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.thermal_solver.main_model_part,self.thermal_solver._get_element_condition_replace_settings()).Execute()

        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

        self.thermal_solver.PrepareModelPart()
        
        for elem in self.thermal_solver.main_model_part.Elements:
            print(f"Element ID: {elem.Id}, Thermal Model Part Elements: {elem.Info()}")
            print(f"Element ID: {elem.Id}, Thermal Model Part Nodes: {elem.GetNodes()}")
        print("Number of Elements in Structural Solver:", len(self.thermal_solver.main_model_part.Elements))
        print("Number of Nodes in Structural Solver:", len(self.thermal_solver.main_model_part.Nodes))
        print("Number of Elements in Structural Solver:", len(self.structural_solver.main_model_part.Elements))
        print("Number of Nodes in Structural Solver:", len(self.structural_solver.main_model_part.Nodes))



    def AddDofs(self):
        self.structural_solver.AddDofs()
        self.thermal_solver.AddDofs()

        #for node in self.structural_solver.main_model_part.Nodes:
            #print(f"Node ID: {node.Id}, Degrees of Freedom: {node.GetDof(KratosMultiphysics.TEMPERATURE)}")

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.structural_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.structural_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.structural_solver.Initialize()
        self.thermal_solver.Initialize()
        #for elem in self.structural_solver.main_model_part.Elements:
            #print(f"Element ID: {elem.Id}, Structural Model Part Elements: {elem.Info()}")
            #print(f"Element ID: {elem.Id}, Structural Model Part Nodes: {elem.GetNodes()}")

    def Clear(self):
        (self.structural_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.structural_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.structural_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

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

        self.structural_solver.InitializeSolutionStep()
        self.structural_solver.Predict()

        KratosMultiphysics.Logger.PrintInfo("\t" + "Solving STRUCTURAL part...")
        structural_is_converged = self.structural_solver.SolveSolutionStep()

        self.RemoveConvectiveVelocity()

        return structural_is_converged and thermal_is_converged

    def FinalizeSolutionStep(self):
        self.structural_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def RemoveConvectiveVelocity(self):
        if self.is_dynamic:
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VELOCITY, KratosMultiphysics.MESH_VELOCITY, self.thermal_solver.GetComputingModelPart(), self.thermal_solver.GetComputingModelPart(), 0)

