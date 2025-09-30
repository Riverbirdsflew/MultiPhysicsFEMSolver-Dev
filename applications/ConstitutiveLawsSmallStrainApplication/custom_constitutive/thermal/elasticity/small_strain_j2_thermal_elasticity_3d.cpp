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
#include "thermal_elastic_isotropic_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/advanced_constitutive_law_small_strain_utilities.h"
#include "includes/checks.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mThermalStrain = ZeroVector(this->GetStrainSize());
    mThermalStrainNormal = 0.0;
}

//************************************************************************************
//************************************************************************************

void ThermalElasticIsotropic3D::FinalizeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    Vector thermal_strain;
    this->CalculateMaterialResponseCauchy(rValues, thermal_strain);
    mThermalStrain = thermal_strain;
    mThermalStrainNormal = thermal_strain[0];
    this->SetThermalStrain(thermal_strain);
    this->SetThermalStrainNormal(thermal_strain[0]);
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    Vector thermal_strain;

    Flags& r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType& r_strain_vector = rValues.GetStrainVector();

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We add the initial strains
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

    // Calculate thermal strain
    thermal_strain.resize(6, false);
    thermal_strain = this->GetThermalStrain();
    Vector ThermalStrainIncerementVector(6);
    this->CalculateThermalStrainIncrement(ThermalStrainIncerementVector, rValues);
    noalias(thermal_strain) += ThermalStrainIncerementVector;

    //substract themal strain
    noalias(r_strain_vector) -= thermal_strain;


    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
        AddInitialStressVectorContribution<StressVectorType>(r_stress_vector);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculateThermalStrainIncrement(
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
    rThermalStrainIncrementVector[4] = 0.0;
    rThermalStrainIncrementVector[5] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

int ThermalElasticIsotropic3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_ERROR_IF_NOT(rElementGeometry[0].SolutionStepsDataHas(TEMPERATURE))  << "The TEMPERATURE variable is not available at the nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT)) << "The THERMAL_EXPANSION_COEFFICIENT is not set in the material properties." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[THERMAL_EXPANSION_COEFFICIENT] < 0.0)   << "The THERMAL_EXPANSION_COEFFICIENT is negative..." << std::endl;
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

double& ThermalElasticIsotropic3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable, double& rValue
    )
{
    if (rThisVariable == REFERENCE_TEMPERATURE) {
        rValue = mReferenceTemperature;
    } else {
        BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
    return (rValue);
}


/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculatePK2Stress(
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

    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculateElasticMatrix(
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
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(rConstitutiveMatrix, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
