# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosConcentrationDiffusionApplication import *

application = KratosConcentrationDiffusionApplication()
application_name = "KratosConcentrationDiffusionApplication"

_ImportApplication(application, application_name)
