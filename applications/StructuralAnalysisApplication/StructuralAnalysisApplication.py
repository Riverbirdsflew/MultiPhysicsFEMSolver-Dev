# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosStructuralAnalysisApplication import *

application = KratosStructuralAnalysisApplication()
application_name = "KratosStructuralAnalysisApplication"

_ImportApplication(application, application_name)

if CheckIfApplicationsAvailable("ConstitutiveLawsSmallStrainApplication"):
    # if available import the advanced constitutive laws for small strain
    import KratosMultiphysics.ConstitutiveLawsSmallStrainApplication

from . import python_registry_lists
python_registry_utilities.RegisterModelersList("KratosMultiphysics.StructuralAnalysisApplication", python_registry_lists)
python_registry_utilities.RegisterOperationsList("KratosMultiphysics.StructuralAnalysisApplication", python_registry_lists)
python_registry_utilities.RegisterProcessesList("KratosMultiphysics.StructuralAnalysisApplication", python_registry_lists)
