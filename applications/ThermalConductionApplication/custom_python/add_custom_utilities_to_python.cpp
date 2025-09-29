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
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/face_heat_utilities.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos::Python {

void GenerateModelPart(FaceHeatUtilities& FaceHeatUtilities,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
{
    if(domain_size == 2)
    {
        FaceHeatUtilities.GenerateModelPart(origin_model_part, destination_model_part, KratosComponents<Element>::Get("ThermCond2D"),KratosComponents<Condition>::Get("ThermalFace2D2N")	);
    }
    else if(domain_size == 3)
    {
        FaceHeatUtilities.GenerateModelPart(origin_model_part, destination_model_part,KratosComponents<Element>::Get("ThermCond3D"),KratosComponents<Condition>::Get("ThermalFace3D3N")	);
    }
}

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<FaceHeatUtilities>(m,"FaceHeatUtilities").def(py::init<>())
    .def("ApplyFaceHeat",&FaceHeatUtilities::ApplyFaceHeat)
    .def("ConditionModelPart",&FaceHeatUtilities::ConditionModelPart)
    .def("GenerateModelPart",GenerateModelPart)
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

}

} // namespace Kratos::Python
