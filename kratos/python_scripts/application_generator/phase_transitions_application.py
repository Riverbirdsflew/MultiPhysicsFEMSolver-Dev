from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

# Set the application name and generate Camel, Caps and Low
appNameCamel = "PhaseTransitions"

# Fetch the applications directory
debugApp = ApplicationGenerator(appNameCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='TRANSFORMATION_STRAIN', vtype='Vector'),
    VariableCreator(name='TRANSFORMATION_STRESS', vtype='Vector'),
    VariableCreator(name='TRANSFORMATION_PLASTICITY_STRAIN', vtype='Vector'),
    VariableCreator(name='TRANSFORMATION_PLASTICITY_STRESS', vtype='Vector'),
    VariableCreator(name='CREEP_STRAIN', vtype='Vector'),
    VariableCreator(name='CREEP_STRESS', vtype='Vector'),
])


# debugApp.AddProcesses([
#     ProcessCreator('DoSomethingProcess')
# ])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appNameCamel))
