//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "therm_diff_struct_phase_application.h"
#include "therm_diff_struct_phase_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos::Python {

PYBIND11_MODULE(KratosThermDiffStructPhaseApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosThermDiffStructPhaseApplication,
        KratosThermDiffStructPhaseApplication::Pointer,
        KratosApplication>(m, "KratosThermDiffStructPhaseApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUSTENITE_MASS_FRACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FERRITE_MASS_FRACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CEMENTITE_MASS_FRACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PEARLITE_MASS_FRACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPPER_BAINITE_MASS_FRACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOWER_BAINITE_MASS_FRACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MARTENSITE_MASS_FRACTION )

    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Kratos::Python
