
import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input should be a Kratos Model object")

    if (type(solver_settings) != KratosMultiphysics.Parameters):
        raise Exception("input should be a Kratos Parameters object")

    solver_type = solver_settings["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "therm_phase_coupled" or solver_type == "ThermPhaseCoupled"):
            solver_module_name = "therm_phase_coupled_solver"
        elif (solver_type == "therm_diff_phase_coupled" or solver_type == "ThermDiffPhaseCoupled"):
            solver_module_name = "therm_diff_phase_coupled_solver"
        elif (solver_type == "therm_struct_phase_coupled" or solver_type == "ThermStructPhaseCoupled"):
            solver_module_name = "therm_struct_phase_coupled_solver"
        elif (solver_type == "therm_diff_struct_phase_coupled" or solver_type == "ThermDiffStructPhaseCoupled"):
            solver_module_name = "therm_diff_struct_phase_coupled_solver"
        else:
            err_msg = 'Requested solver_type: ' + solver_type + ' is not available.'
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "therm_phase_coupled" or solver_type == "ThermPhaseCoupled"):
            solver_module_name = "trilinos_therm_phase_coupled_solver"
        elif (solver_type == "therm_diff_phase_coupled" or solver_type == "ThermDiffPhaseCoupled"):
            solver_module_name = "trilinos_therm_diff_phase_coupled_solver"
        elif (solver_type == "therm_struct_phase_coupled" or solver_type == "ThermStructPhaseCoupled"):
            solver_module_name = "trilinos_therm_struct_phase_coupled_solver"
        elif (solver_type == "therm_diff_struct_phase_coupled" or solver_type == "ThermDiffStructPhaseCoupled"):
            solver_module_name = "trilinos_therm_diff_struct_phase_coupled_solver"
        else:
            err_msg = 'Requested solver_type: ' + solver_type + ' is not available.'
            raise Exception(err_msg)
    
    # Wrong parallelism check
    else:
        err_msg = "Parallelism is neither OpenMP nor MPI."
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.ThermDiffStructPhaseApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input should be a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input should be a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_settings = custom_settings["solver_settings"]

    return CreateSolverByParameters(model, solver_settings, parallelism)