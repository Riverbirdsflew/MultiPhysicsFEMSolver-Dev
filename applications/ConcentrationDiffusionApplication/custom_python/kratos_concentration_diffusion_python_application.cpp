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

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "concentration_diffusion_application.h"
#include "concentration_diffusion_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
//#include "custom_python/add_custom_response_functions_to_python.h"//不加入responsefunction


namespace Kratos::Python {

PYBIND11_MODULE(KratosConcentrationDiffusionApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosConcentrationDiffusionApplication,
        KratosConcentrationDiffusionApplication::Pointer,
        KratosApplication>(m, "KratosConcentrationDiffusionApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    //AddCustomResponseFunctionsToPython(m);

    //registering variables in python
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POINT_CONCENTRATION_SOURCE )
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AMBIENT_CONCENTRATION )
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FACE_CONCENTRATION_FLUX )
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONCENTRATION_FLUX )
	  //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONCENTRATION_FACE )
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PHASE_TRANSITION_FRACTION )
	  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PHASE_NAME )

    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Kratos::Python

#endif // KRATOS_PYTHON defined
