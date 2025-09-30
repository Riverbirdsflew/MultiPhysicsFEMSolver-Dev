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
#include "small_strain_j2_thermal_elasticity_plane_stress.h"
#include "includes/checks.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "therm_diff_struct_phase_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/advanced_constitutive_law_small_strain_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStress::CalculatePK2Stress(
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
    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStress(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStress::CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues, Matrix &rElasticityTensor)
{
    
    const double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    const double poisson_ratio = rMaterialProperties.GetValue(POISSON_RATIO);
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticityTensor.size1() != 4 || rElasticityTensor.size2() != 4)
        rElasticityTensor.resize(4, 4, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2. * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = 0.;
    rElasticityTensor(0, 3) = 0.;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2. * mu;
    rElasticityTensor(1, 2) = 0.;
    rElasticityTensor(1, 3) = 0.;
    rElasticityTensor(2, 0) = 0.;
    rElasticityTensor(2, 1) = 0.;
    rElasticityTensor(2, 2) = 0.;
    rElasticityTensor(2, 3) = 0.;
    rElasticityTensor(3, 0) = 0.;
    rElasticityTensor(3, 1) = 0.;
    rElasticityTensor(3, 2) = 0.;
    rElasticityTensor(3, 3) = mu;
}


/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStress::CalculateThermalStrainIncrement(
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

    rThermalStrainIncrementVector[0] = alpha * delta_temperature;
    rThermalStrainIncrementVector[1] = alpha * delta_temperature;
    rThermalStrainIncrementVector[2] = alpha * delta_temperature;
    rThermalStrainIncrementVector[3] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ThermalElasticityPlaneStress::CalculateCauchyGreenStrain(
    Parameters& rValues,
    ConstitutiveLaw::StrainVectorType& rStrainVector
    )
{
    ConstitutiveLawUtilities<3>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
