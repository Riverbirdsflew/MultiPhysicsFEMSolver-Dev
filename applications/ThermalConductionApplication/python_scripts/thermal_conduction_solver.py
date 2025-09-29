# Importing the Kratos Library
import KratosMultiphysics

# Auxiliary function to check the parallel type at runtime
#TODO: Delete this once we come up with the final factory-based design
def _CheckIsDistributed():
    if KratosMultiphysics.ParallelEnvironment.HasDataCommunicator("World"):
        world_data_comm = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")
        return world_data_comm.IsDistributed()
    else:
        return False

# Import applications
import KratosMultiphysics.ThermalConductionApplication

# If required, import parallel applications and modules
if _CheckIsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    import KratosMultiphysics.mpi.distributed_import_model_part_utility as distributed_import_model_part_utility

# Importing factories
if _CheckIsDistributed():
    import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as linear_solver_factory
else:
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    import KratosMultiphysics.base_convergence_criteria_factory as convergence_criteria_factory

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics import auxiliary_solver_utilities

def CreateSolver(model, custom_settings):
    return ThermalConductionSolver(model, custom_settings)
    

class ThermalConductionSolver(PythonSolver):
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

        KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "Construction finished")

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

        if target_model_part == None:
            target_model_part = self.main_model_part

        ''' Add nodal solution step variables
        '''
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LATENT_AUSTENIZE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LATENT_EUTECTOID_DECOMPOSITION)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LATENT_BAINITE_TRANSFORMATION)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LATENT_MARTENSITE_TRANSFORMATION)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_AUSTENITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_FERRITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_CEMENTITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_PEARLITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_UPPER_BAINITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_LOWER_BAINITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MASS_FRACTION_MARTENSITE)
        target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_FLUX)


        # Adding nodal area variable (some solvers use it. TODO: Ask)
        #target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        # If LaplacianElement is used
        if (self.settings["element_replace_settings"]["element_name"].GetString() == "LaplacianElement"):
            target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        # If MPI distributed, add the PARTITION_INDEX
        if _CheckIsDistributed():
            target_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        auxiliary_solver_utilities.AddVariables(self.main_model_part, self.settings["auxiliary_variables_list"])

        KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def AddDofs(self):
        # Set DOFs and reaction variables list from Kratos parameters settings
        dofs_with_reactions_list = []
        dof_var_name = "TEMPERATURE"
        reaction_var_name = "REACTION_FLUX"
        dofs_with_reactions_list.append([dof_var_name,reaction_var_name])
        if self.settings["gradient_dofs"].GetBool():
            grad_dof_var_name = "TEMPERATURE_GRADIENT"
            grad_react_var_name = "REACTION"
            comp_list = ["_X","_Y"] if self.settings["domain_size"].GetInt() == 2 else ["_X","_Y","_Z"]
            for comp in comp_list:
                dofs_with_reactions_list.append([grad_dof_var_name+comp,grad_react_var_name+comp])

        # Add the DOFs and reaction list to each node
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_with_reactions_list, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "DOF's ADDED")

    def GetDofsList(self):
        """This function creates and returns a list with the DOFs defined in the Kratos parameters settings
        Note that element GetSpecifications method cannot be used in this case as DOF variables are a priori unknown
        """
        dofs_list = []
        dofs_list.append("TEMPERATURE")
        if self.settings["gradient_dofs"].GetBool():
            grad_dof_var_name = "TEMPERATURE_GRADIENT"
            comp_list = ["_X","_Y"] if self.settings["domain_size"].GetInt() == 2 else ["_X","_Y","_Z"]
            for comp in comp_list:
                dofs_list.append(grad_dof_var_name + comp)

        return dofs_list

    def ImportModelPart(self):
        """This function imports the ModelPart"""
        if self.solver_imports_model_part:
            if not _CheckIsDistributed():
                self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])
            else:
                self.distributed_model_part_importer = distributed_import_model_part_utility.DistributedImportModelPartUtility(
                    self.main_model_part,
                    self.settings)
                self.distributed_model_part_importer.ImportModelPart()

    def PrepareModelPart(self):
        assign_neighbour_elements = self.settings["assign_neighbour_elements_to_conditions"].GetBool()
        if not self.is_restarted():
            # Import material properties
            materials_imported = self.import_materials()
            if materials_imported:
                KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "Materials were successfully imported.")
            else:
                KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "Materials were not imported.")

            # KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part,self._get_element_condition_replace_settings()).Execute()

            # tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
            # throw_errors = False
            # flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
            # if assign_neighbour_elements:
            #     flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
            # else:
            #     flags |= (tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS).AsFalse()
            # tmoc(self.main_model_part,throw_errors, flags).Execute()

            self._set_and_fill_buffer()

        # Create the MPI communicators
        if _CheckIsDistributed():
            self.distributed_model_part_importer.CreateCommunicators()

        if (self.settings["echo_level"].GetInt() > 0):
            KratosMultiphysics.Logger.PrintInfo(self.model)

        KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part."""
        KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "Initializing ...")
        # The thermal_conduction solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        thermal_conduction_solution_strategy = self._GetSolutionStrategy()
        thermal_conduction_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        if not self.is_restarted():
            thermal_conduction_solution_strategy.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of ImplicitSolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            try:
                thermal_conduction_solution_strategy.SetInitializePerformedFlag(True)
            except AttributeError:
                pass
        self.Check()
        KratosMultiphysics.Logger.PrintInfo("::[ThermalConductionSolver]:: ", "Finished initialization.")

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
        if not self.main_model_part.IsDistributed():
            convergence_criterion = convergence_criteria_factory.ConvergenceCriteriaFactory(self._get_convergence_criterion_settings())
            return convergence_criterion.convergence_criterion
        else:
            convergence_criterion = self.__base_convergence_criteria_factory_mpi(self._get_convergence_criterion_settings())
            return convergence_criterion

    def _CreateLinearSolver(self):
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        if not self.main_model_part.IsDistributed():
            # Set the serial builder and solver
            if self.settings["block_builder"].GetBool():
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        else:
            # Set the Epetra vectors communicator
            epetra_communicator = self.get_epetra_communicator()

            # Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
            domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            if domain_size == 3:
                guess_row_size = 20
            else:
                guess_row_size = 10

            # Set the parallel builder and solver
            if self.settings["block_builder"].GetBool():
                builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolver(
                    epetra_communicator,
                    guess_row_size,
                    linear_solver)
            else:
                builder_and_solver = KratosTrilinos.TrilinosEliminationBuilderAndSolver(
                    epetra_communicator,
                    guess_row_size,
                    linear_solver)

        return builder_and_solver

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
        computing_model_part = self.GetComputingModelPart()
        thermal_conduction_scheme = self._GetScheme()
        thermal_conduction_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        if not computing_model_part.IsDistributed():
            return KratosMultiphysics.LineSearchStrategy(
                computing_model_part,
                thermal_conduction_scheme,
                thermal_conduction_convergence_criterion,
                builder_and_solver,
                self.settings["max_iteration"].GetInt(),
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["move_mesh_flag"].GetBool())
        else:
            err_msg = "\'line_search\' solution strategy is not MPI compatible."
            raise Exception(err_msg)

    #TODO: THIS MUST BE IMPLEMENTED IN A base_convergence_criteria_factory_mpi.py
    #TODO: THEN WE CAN IMPORT IT AS convergence_criteria_factory TO AVOID DISTINGUISHING THE SERIAL AND THE PARALLEL FACTORIES
    def __base_convergence_criteria_factory_mpi(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        D_RT = convergence_criterion_parameters["solution_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["solution_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()

        echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        convergence_crit = convergence_criterion_parameters["convergence_criterion"].GetString()

        if(echo_level >= 1):
            KratosMultiphysics.Logger.PrintInfo("::[ConvergenceCriterionFactory]:: ", "CONVERGENCE CRITERION : " +
                convergence_criterion_parameters["convergence_criterion"].GetString())

        if(convergence_crit == "solution_criterion"):
            convergence_criterion = KratosTrilinos.TrilinosDisplacementCriteria(D_RT, D_AT)
            convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "residual_criterion"):
            convergence_criterion = KratosTrilinos.TrilinosResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(echo_level)

        elif(convergence_crit == "and_criterion"):
            Displacement = KratosTrilinos.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosTrilinos.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosTrilinos.TrilinosAndCriteria(Residual, Displacement)

        elif(convergence_crit == "or_criterion"):
            Displacement = KratosTrilinos.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosTrilinos.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosTrilinos.TrilinosOrCriteria(Residual, Displacement)
        else:
            err_msg =  "The requested convergence criterion \"" + convergence_crit + "\" is not available!\n"
            err_msg += "Available options are: \"solution_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
            raise Exception(err_msg)

        return convergence_criterion

    def get_epetra_communicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = KratosTrilinos.CreateCommunicator()
        return self._epetra_communicator

