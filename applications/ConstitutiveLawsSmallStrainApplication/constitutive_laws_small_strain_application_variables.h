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
//  Main authors:    whf
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/mat_variables.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{

    enum class TangentOperatorEstimation {Analytic = 0, FirstOrderPerturbation = 1, SecondOrderPerturbation = 2, Secant = 3, SecondOrderPerturbationV2 = 4, InitialStiffness = 5, OrthogonalSecant = 6};

///@name  Functions
///@{
    // Constitutive laws variables

    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, int, MAX_NUMBER_NL_CL_ITERATIONS)

    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, double, EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, Matrix, PLASTIC_STRAIN_TENSOR)
    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, Vector, PLASTIC_STRAIN_VECTOR)

    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, double, VON_MISES_STRESS)
    
    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, double, YIELD_STRENGTH)

    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, Vector, THERMAL_STRAIN_VECTOR)
    KRATOS_DEFINE_APPLICATION_VARIABLE(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION, double, THERMAL_STRAIN_NORMAL)

}
