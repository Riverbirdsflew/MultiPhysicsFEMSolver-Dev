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
#include "concentration_diffusion_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( double, POINT_CONCENTRATION_SOURCE )
KRATOS_CREATE_VARIABLE( double, AMBIENT_CONCENTRATION )
KRATOS_CREATE_VARIABLE( double, FACE_CONCENTRATION_FLUX )
KRATOS_CREATE_VARIABLE( double, CONCENTRATION_FLUX )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONCENTRATION_GRADIENT)

KRATOS_CREATE_VARIABLE( double, PHASE_TRANSITION_FRACTION )
KRATOS_CREATE_VARIABLE( double, PHASE_NAME )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

KRATOS_CREATE_VARIABLE(double, PROJECTED_SCALAR1)
KRATOS_CREATE_VARIABLE(double, TRANSFER_COEFFICIENT)

}
