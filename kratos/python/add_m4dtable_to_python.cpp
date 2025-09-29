//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marx Xu
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/m4dtable.h"
#include "add_m4dtable_to_python.h"

namespace Kratos::Python
{

namespace py = pybind11;

typedef Kratos::M4dTable M4dTableType;

void M4dTableGetNearestValue(M4dTableType& ThisM4dTable, Kratos::M4dTable::RowType const& m, double const& x, double& lx, double& rx)
{
    Kratos::M4dTable::ColumnType ln, rn;

    if(m < ThisM4dTable.Dim())
    {
        Kratos::M4dTable::ValueContainerType& arg = ThisM4dTable.Argument()[m].second;

        if(arg.size() > 0)
        {
            ThisM4dTable.GetNearestColumn(m, x, ln, rn);
            lx = arg[ln];
            rx = arg[rn];
        }
    }
}

void  AddM4dTableToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<M4dTableType, M4dTableType::Pointer>(m,"Multi4DTable")
    .def(py::init<>())
    .def(py::init<Kratos::M4dTable::RowType const&>())
    .def("GetValue", &M4dTableType::GetInterpolateValue)
    .def("GetDerivative",&M4dTableType::GetDerivativeValue)
    .def("GetIntegration",&M4dTableType::GetIntegrationValue)
    .def("GetNearestValue", M4dTableGetNearestValue)
    .def("Dim", static_cast<Kratos::M4dTable::RowType& (M4dTableType::*)()>(&M4dTableType::Dim))
    .def("AddX", static_cast<void (M4dTableType::*)(Kratos::M4dTable::RowType const&, Kratos::M4dTable::ValueType const&)>(&M4dTableType::PushBack))
    .def("AddY", static_cast<void (M4dTableType::*)(Kratos::M4dTable::ValueType const&)>(&M4dTableType::PushBack))
    .def("Clear", &M4dTableType::Clear)
    .def("__str__", PrintObject<M4dTableType>)
    ;
}

}  // namespace Kratos::Python.

