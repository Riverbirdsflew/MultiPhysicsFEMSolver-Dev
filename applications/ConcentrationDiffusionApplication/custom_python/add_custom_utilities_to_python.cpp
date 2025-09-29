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
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/face_concentration_utilities.h"
//#include "custom_utilities/pure_concentration_convection_tools.h"
//#include "custom_utilities/pure_concentration_convection_CrankN_tools.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{
	
void GenerateModelPart(FaceConcentrationUtilities& FaceConcentrationUtilities,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
{
    if(domain_size == 2)
    {
        FaceConcentrationUtilities.GenerateModelPart(origin_model_part, destination_model_part, KratosComponents<Element>::Get("ConcenDiff2D"),KratosComponents<Condition>::Get("ConcentrationFace2D2N")	);
    }
    else if(domain_size == 3)
    {
        FaceConcentrationUtilities.GenerateModelPart(origin_model_part, destination_model_part,KratosComponents<Element>::Get("ConcenDiff3D"),KratosComponents<Condition>::Get("ConcentrationFace3D3N")	);
    }
}



  void  AddCustomUtilitiesToPython(pybind11::module& m)
  {
    namespace py = pybind11;
    
    py::class_<FaceConcentrationUtilities>(m,"FaceConcentrationUtilities").def(py::init<>())
    .def("ApplyFaceConcentration",&FaceConcentrationUtilities::ApplyFaceConcentration)
    .def("ConditionModelPart",&FaceConcentrationUtilities::ConditionModelPart)
    .def("GenerateModelPart",GenerateModelPart)
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    //~ py::class_< PureConcentrationConvectionUtilities< 2, SparseSpaceType, LinearSolverType >>(m,"PureConcentrationConvectionUtilities2D").def(py::init<	>() )
    //~ .def("ConstructSystem",&PureConcentrationConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
    //~ .def("CalculateProjection",&PureConcentrationConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::CalculateProjection)
    //~ .def("ConvectScalarVar",&PureConcentrationConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    //~ .def("ClearSystem",&PureConcentrationConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
    //~ ;

    //~ py::class_< PureConcentrationConvectionUtilities< 3, SparseSpaceType, LinearSolverType >>(m,"PureConcentrationConvectionUtilities3D").def(py::init<	>() )
    //~ .def("ConstructSystem",&PureConcentrationConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
    //~ .def("CalculateProjection",&PureConcentrationConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::CalculateProjection)
    //~ .def("ConvectScalarVar",&PureConcentrationConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    //~ .def("ClearSystem",&PureConcentrationConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
    //~ ;

    //~ py::class_< PureConcentrationConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >>(m,"PureConcentrationConvectionCrankNUtilities2D").def(py::init<	>() )
    //~ .def("ConstructSystem",&PureConcentrationConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
    //~ .def("CalculateProjection",&PureConcentrationConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::CalculateProjection)
    //~ .def("ConvectScalarVar",&PureConcentrationConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    //~ .def("ClearSystem",&PureConcentrationConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
    //~ ;

    //~ py::class_< PureConcentrationConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >>(m,"PureConcentrationConvectionCrankNUtilities3D").def(py::init<	>() )
    //~ .def("ConstructSystem",&PureConcentrationConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
    //~ .def("CalculateProjection",&PureConcentrationConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::CalculateProjection)
    //~ .def("ConvectScalarVar",&PureConcentrationConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    //~ .def("ClearSystem",&PureConcentrationConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
    //~ ;


  }

}  // namespace Python.

} // Namespace Kratos
