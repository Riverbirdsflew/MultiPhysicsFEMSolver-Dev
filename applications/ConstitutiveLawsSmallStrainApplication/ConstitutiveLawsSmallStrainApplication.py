# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.StructuralAnalysisApplication # temp until dependencies are moved to Core
from KratosConstitutiveLawsSmallStrainApplication import *
application = KratosConstitutiveLawsSmallStrainApplication()
application_name = "KratosConstitutiveLawsSmallStrainApplication"

_ImportApplication(application, application_name)
