
import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.ThermDiffStructPhaseApplication import python_solvers_wrapper_therm_diff_struct_phase

class ThermDiffStructPhaseAnalysis(AnalysisStage):
    '''Main script for ThermDiffStructPhase simulations.'''

    def __init__(self, model, project_parameters):
        
        super(ThermDiffStructPhaseAnalysis, self).__init__(model, project_parameters)

    def _CreateSolver(self):
        parallelism = self.project_parameters["problem_data"]["parallel_type"].GetString()
        return python_solvers_wrapper_therm_diff_struct_phase.CreateSolver(self.model, self.project_parameters, parallelism)

    def _GetSimulationName(self):
        if not self.project_parameters["problem_data"].Has("problem_name"):
            return "::[Therm-Diff-Struct-Phase Simulation]:: "
        else:
            return self.project_parameters["problem_data"]["problem_name"].GetString()

    def _GetOrderOfProcessesInitialization(self):
        return ["constraints_process_list","loads_process_list", "initial_conditions_process_list"]

    def _GetOrderOfOutputProcessesInitialization(self):
        return ["gid_output","vtk_output"]

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fsi_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fsi_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ThermDiffStructPhaseAnalysis(model, parameters)
    simulation.Run()
