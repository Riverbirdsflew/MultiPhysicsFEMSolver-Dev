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
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
//#include "custom_strategies/strategies/residualbased_eulerian_concendiff_strategy.h"
//#include "custom_strategies/strategies/residualbased_semi_eulerian_concendiff_strategy.h"
//#include "custom_strategies/strategies/explicit_runge_kutta_4_eulerian_concendiff_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//convergence criterias
#include "custom_strategies/strategies/residualbased_concendiff_strategy.h"
#include "custom_strategies/strategies/residualbased_concendiff_strategy_nonlinear.h"


namespace Kratos::Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    //typedef ExplicitSolvingStrategyRungeKutta4< SparseSpaceType, LocalSpaceType > ExplicitSolvingStrategyRungeKutta4Type;
    //typedef ExplicitBuilder< SparseSpaceType, LocalSpaceType > ExplicitBuilderType;

    
    py::class_< ResidualBasedConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            ResidualBasedConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"ResidualBasedConcentrationDiffusionStrategy")
            .def(py::init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int	>() )
            .def("Clear",&ResidualBasedConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;

    //~ py::class_< ResidualBasedEulerianConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            //~ ResidualBasedEulerianConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            //~ BaseSolvingStrategyType >
            //~ (m,"ResidualBasedEulerianConcentrationDiffusionStrategy")
            //~ .def(py::init<	ModelPart&, LinearSolverType::Pointer,	bool, int	>() )
            //~ .def("Clear",&ResidualBasedEulerianConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            //~ ;

    //~ py::class_< ResidualBasedSemiEulerianConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            //~ ResidualBasedSemiEulerianConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            //~ BaseSolvingStrategyType >
            //~ (m,"ResidualBasedSemiEulerianConcentrationDiffusionStrategy")
            //~ .def(py::init<	ModelPart&, LinearSolverType::Pointer,	bool, int	>() )
            //~ .def("Clear",&ResidualBasedSemiEulerianConcentrationDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            //~ ;

    py::class_< ResidualBasedConcentrationDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            ResidualBasedConcentrationDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"ResidualBasedConcentrationDiffusionStrategyNonLinear")
            .def(py::init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int ,double	>() )
            .def("Clear",&ResidualBasedConcentrationDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;

    //~ typedef ExplicitSolvingStrategyRungeKutta4ConcentrationDiffusion< SparseSpaceType, LocalSpaceType > ExplicitSolvingStrategyRungeKutta4ConcentrationDiffusionType;
    //~ py::class_<ExplicitSolvingStrategyRungeKutta4ConcentrationDiffusionType, typename ExplicitSolvingStrategyRungeKutta4ConcentrationDiffusionType::Pointer, ExplicitSolvingStrategyRungeKutta4Type>(m, "ExplicitSolvingStrategyRungeKutta4ConcentrationDiffusion")
        //~ .def(py::init<ModelPart &, bool, int>())
        //~ .def(py::init<ModelPart&, typename ExplicitBuilderType::Pointer, bool, int>())
        //~ ;

}

} // namespace Kratos::Python
