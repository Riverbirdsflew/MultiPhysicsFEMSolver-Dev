# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ThermalConductionApplication

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics import auxiliary_solver_utilities

def CreateSolver(model, custom_settings):
    return PhaseTransitionSolver(model, custom_settings)
    

class PhaseTransitionSolver(PythonSolver):
    """The base class for thermal-conduction solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _CreateScheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _CreateScheme
    _CreateConvergenceCriterion
    _CreateLinearSolver
    _CreateBuilderAndSolver
    _CreateSolutionStrategy

    The thermal_conduction_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions _GetSolutionStrategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

        # Thermal Conduction variables check
        #self._ThermalConductionVariablesCheck(custom_settings)

        model_part_name = self.settings["model_part_name"].GetString()

        # Set default buffer size
        self.min_buffer_size = 1

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
            self.solver_imports_model_part = False
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name) # Model.CreateodelPart()
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.solver_imports_model_part = True

        KratosMultiphysics.Logger.PrintInfo("::[PhaseTransitionSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "ThermalModelPart",
            "domain_size" : -1,
            "echo_level": 0,
            "analysis_type": "linear",
            "solver_type": "thermal_conduction_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping" : {
                "time_step": 1.0
            },
            "calculation_flag_settings"		:{
                "is_phase_latent_calculated" : false,
                "is_thermal_strain_calculated": false,
                "is_phase_transition_strain_calculated" : false,
                "is_trip_calculated" : false,
                "is_creep_calculated" : false,
                "is_thermo_calculated" : false,
                "is_structure_calculated" : false,
                "is_phase_transition_calculated" : false,
                "is_diffusion_calculated" : false
            },
            "reform_dofs_at_each_step": false,
            "gradient_dofs" : false,
            "line_search": false,
            "compute_reactions": true,
            "block_builder": true,
            "clear_storage": false,
            "move_mesh_flag": false,
            "convergence_criterion": "residual_criterion",
            "solution_relative_tolerance": 1.0e-4,
            "solution_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "amgcl",
                "smoother_type":"ilu0",
                "krylov_type":"gmres",
                "coarsening_type":"aggregation",
                "max_iteration": 5000,
                "tolerance": 1e-9,
                "scaling": false
            },
            "element_replace_settings" : {
                "element_name" : "EulerianThermCondElement",
                "condition_name" : "ThermalFace"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "auxiliary_variables_list" : [],
            "assign_neighbour_elements_to_conditions" : true
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def AddVariables(self, target_model_part=None):

        pass

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def AddDofs(self):
        pass

    def GetDofsList(self):
        pass

    def ImportModelPart(self):
        pass

    def PrepareModelPart(self):
        pass

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part."""
        KratosMultiphysics.Logger.PrintInfo("::[PhaseTransitionSolver]:: ", "Initializing ...")
        # The thermal_conduction solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        KratosMultiphysics.IS_PHASE_LATENT_CALCULATED = self.settings["calculation_flag_settings"]["is_phase_latent_calculated"].GetBool()
        KratosMultiphysics.IS_THERMAL_STRAIN_CALCULATED = self.settings["calculation_flag_settings"]["is_thermal_strain_calculated"].GetBool()
        KratosMultiphysics.IS_PHASE_TRANSITION_STRAIN_CALCULATED = self.settings["calculation_flag_settings"]["is_phase_transition_strain_calculated"].GetBool()
        KratosMultiphysics.IS_TRIP_CALCULATED = self.settings["calculation_flag_settings"]["is_trip_calculated"].GetBool()
        KratosMultiphysics.IS_CREEP_CALCULATED = self.settings["calculation_flag_settings"]["is_creep_calculated"].GetBool()
        KratosMultiphysics.IS_THERMO_CALCULATED = self.settings["calculation_flag_settings"]["is_thermo_calculated"].GetBool()
        KratosMultiphysics.IS_STRUCTURE_CALCULATED = self.settings["calculation_flag_settings"]["is_structure_calculated"].GetBool()
        KratosMultiphysics.IS_PHASE_TRANSITION_CALCULATED = self.settings["calculation_flag_settings"]["is_phase_transition_calculated"].GetBool()
        KratosMultiphysics.IS_DIFFUSION_CALCULATED = self.settings["calculation_flag_settings"]["is_diffusion_calculated"].GetBool()

        KratosMultiphysics.Logger.PrintInfo("::[PhaseTransitionSolver]:: ", "Finished initialization.")

    def GetOutputVariables(self):
        pass

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        thermal_conduction_solution_strategy = self._GetSolutionStrategy()
        thermal_conduction_solution_strategy.Solve()

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def GetComputingModelPart(self):
        return self.main_model_part

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self._GetSolutionStrategy().SetEchoLevel(level)

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def Check(self):
        self._GetSolutionStrategy().Check()

    #### Specific internal functions ####

    def _GetScheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._CreateScheme()
        return self._solution_scheme

    def _GetConvergenceCriterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._CreateConvergenceCriterion()
        return self._convergence_criterion

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_thermal_conduction_solution_strategy'):
            self._thermal_conduction_solution_strategy = self._CreateSolutionStrategy()
        return self._thermal_conduction_solution_strategy

    def import_materials(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)

            # We set the properties that are nodal
            self._assign_nodally_properties()
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    #### Private functions ####

    def _assign_nodally_properties(self):

        # We transfer the values of the con.diff variables to the nodes
        with open(self.settings["material_import_settings"]["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KratosMultiphysics.Parameters(parameter_file.read())

        for i in range(materials["properties"].size()):
            model_part = self.model.GetModelPart(materials["properties"][i]["model_part_name"].GetString())
            mat = materials["properties"][i]["Material"]

            for key, value in mat["Variables"].items():
                var = KratosMultiphysics.KratosGlobals.GetVariable(key)
                if value.IsDouble():
                    KratosMultiphysics.VariableUtils().SetVariable(var, value.GetDouble(), model_part.Nodes)
                elif value.IsVector():
                    KratosMultiphysics.VariableUtils().SetVariable(var, value.GetVector(), model_part.Nodes)
                else:
                    raise ValueError("Type of value is not available")

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information."""
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        required_buffer_size = self.GetMinimumBufferSize()
        current_buffer_size = self.main_model_part.GetBufferSize()
        buffer_size = max(current_buffer_size, required_buffer_size)
        self.main_model_part.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for _ in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            self.main_model_part.CloneTimeStep(time)

    def _get_element_condition_replace_settings(self):
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
        element_list = ["EulerianThermCondElement","LaplacianElement","AxisyEulerianThermCond"]
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
            KratosMultiphysics.Logger.PrintInfo("Replacing Element: ", self.settings["element_replace_settings"]["element_name"])

        ## Conditions
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()
        condition_list = ["FluxCondition","ThermalFace","AxisyThermalFace","LineCondition","SurfaceCondition"]
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
            KratosMultiphysics.Logger.PrintInfo("Replacing Condition: ", self.settings["element_replace_settings"]["condition_name"])

        return self.settings["element_replace_settings"]

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("solution_relative_tolerance",self.settings["solution_relative_tolerance"])
        conv_params.AddValue("solution_absolute_tolerance",self.settings["solution_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _CreateConvergenceCriterion(self):
        pass

    def _CreateLinearSolver(self):
        pass

    def _CreateBuilderAndSolver(self):
        pass

    @classmethod
    def _CreateScheme(self):
        """Create the solution scheme for the thermal-conduction problem."""
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            thermal_conduction_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                thermal_conduction_solution_strategy = self._create_newton_raphson_strategy()
            else:
                thermal_conduction_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return thermal_conduction_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        thermal_conduction_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        if not computing_model_part.IsDistributed():
            return KratosMultiphysics.ResidualBasedLinearStrategy(
                computing_model_part,
                thermal_conduction_scheme,
                builder_and_solver,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                False,
                self.settings["move_mesh_flag"].GetBool())
        else:
            return KratosTrilinos.TrilinosLinearStrategy(
                computing_model_part,
                thermal_conduction_scheme,
                builder_and_solver,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                False,
                self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        thermal_conduction_scheme = self._GetScheme()
        thermal_conduction_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        if not computing_model_part.IsDistributed():
            return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
                computing_model_part,
                thermal_conduction_scheme,
                thermal_conduction_convergence_criterion,
                builder_and_solver,
                self.settings["max_iteration"].GetInt(),
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["move_mesh_flag"].GetBool())
        else:
            return KratosTrilinos.TrilinosNewtonRaphsonStrategy(
                computing_model_part,
                thermal_conduction_scheme,
                thermal_conduction_convergence_criterion,
                builder_and_solver,
                self.settings["max_iteration"].GetInt(),
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["move_mesh_flag"].GetBool())

    def _create_line_search_strategy(self):
        pass

    def get_epetra_communicator(self):
        pass

