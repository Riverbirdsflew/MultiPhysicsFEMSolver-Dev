// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "constitutive_laws_small_strain_application.h"
#include "constitutive_laws_small_strain_application_variables.h"

#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosConstitutiveLawsSmallStrainApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosConstitutiveLawsSmallStrainApplication,
        KratosConstitutiveLawsSmallStrainApplication::Pointer,
        KratosApplication>(m, "KratosConstitutiveLawsSmallStrainApplication")
        .def(py::init<>())
        ;

    AddCustomConstitutiveLawsToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    // Constitutive laws variables

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MAX_NUMBER_NL_CL_ITERATIONS)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PLASTIC_STRAIN_TENSOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PLASTIC_STRAIN_VECTOR)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, VON_MISES_STRESS)
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, YIELD_STRENGTH)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_STRAIN_VECTOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_STRAIN_NORMAL)

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
