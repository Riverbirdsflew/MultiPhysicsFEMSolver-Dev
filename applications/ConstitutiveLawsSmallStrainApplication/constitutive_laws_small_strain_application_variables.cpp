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

// System includes

// External includes

// Project includes
#include "constitutive_laws_small_strain_application_variables.h"

namespace Kratos
{
    // Constitutive laws variables

    KRATOS_CREATE_VARIABLE( int, MAX_NUMBER_NL_CL_ITERATIONS)

    KRATOS_CREATE_VARIABLE( double, EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_CREATE_VARIABLE( Matrix, PLASTIC_STRAIN_TENSOR)
    KRATOS_CREATE_VARIABLE( Vector, PLASTIC_STRAIN_VECTOR)

    KRATOS_CREATE_VARIABLE( double, VON_MISES_STRESS)
    
    KRATOS_CREATE_VARIABLE( double, YIELD_STRENGTH)

    KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR)
    KRATOS_CREATE_VARIABLE( double, THERMAL_STRAIN_NORMAL)

}
