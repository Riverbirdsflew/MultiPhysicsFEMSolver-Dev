from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

# Set the application name and generate Camel, Caps and Low
appNameCamel = "ThermDiffStructPhaseApplication"

# Fetch the applications directory
debugApp = ApplicationGenerator(appNameCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='AUSTENITE_MASS_FRACTION', vtype='double'),
    VariableCreator(name='FERRITE_MASS_FRACTION', vtype='double'),
    VariableCreator(name='CEMENTITE_MASS_FRACTION', vtype='double'),
    VariableCreator(name='PEARLITE_MASS_FRACTION', vtype='double'),
    VariableCreator(name='UPPER_BAINITE_MASS_FRACTION', vtype='double'),
    VariableCreator(name='LOWER_BAINITE_MASS_FRACTION', vtype='double'),
    VariableCreator(name='MARTENSITE_MASS_FRACTION', vtype='double'),
])


# debugApp.AddProcesses([
#     ProcessCreator('DoSomethingProcess')
# ])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appNameCamel))
