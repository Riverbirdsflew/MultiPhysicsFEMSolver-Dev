// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_analysis_application/license.txt
//
//  Main authors:    whf
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "small_strain_j2_thermal_elasticity_plane_strain.h"
#include "includes/checks.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "therm_diff_struct_phase_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/advanced_constitutive_law_small_strain_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStrain::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();

    double E = mpPhaseTransSystem->GetElasticityModulus(r_properties, r_geom, r_N, r_process_info);
    double NU = mpPhaseTransSystem->GetPoisson(r_properties, r_geom, r_N, r_process_info);
    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStrain::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();

    double E = mpPhaseTransSystem->GetElasticityModulus(r_properties, r_geom, r_N, r_process_info);
    double NU = mpPhaseTransSystem->GetPoisson(r_properties, r_geom, r_N, r_process_info);
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStrain(rConstitutiveMatrix, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStrain::CalculateThermalStrainIncrement(
    Vector& rThermalStrainIncrementVector, 
    ConstitutiveLaw::Parameters &rValues
    )
{
    const double previous_temperature_gp = AdvCLutils::CalculateInGaussPoint(TEMPERATURE, rValues, 1);
    const double current_temperature_gp = AdvCLutils::CalculateInGaussPoint(TEMPERATURE, rValues, 0);
    const double delta_temperature = current_temperature_gp - previous_temperature_gp;

    const Properties& r_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();

    double alpha = mpPhaseTransSystem->GetExpansion(r_properties, r_geom, r_N, r_process_info);
    double NU = mpPhaseTransSystem->GetPoisson(r_properties, r_geom, r_N, r_process_info);
    alpha *= (1.0 + NU );

    rThermalStrainIncrementVector[0] = alpha * delta_temperature;
    rThermalStrainIncrementVector[1] = alpha * delta_temperature;
    rThermalStrainIncrementVector[2] = 0.0;
    rThermalStrainIncrementVector[3] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStrain::CalculateCauchyGreenStrain(
    Parameters& rValues,
    ConstitutiveLaw::StrainVectorType& rStrainVector
    )
{
    ConstitutiveLawUtilities<3>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
