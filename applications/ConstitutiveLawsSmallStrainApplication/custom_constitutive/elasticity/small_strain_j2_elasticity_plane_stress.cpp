// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| ANALYSIS
//
//  License:         BSD License
//                   license: constitutive_laws_small_strain_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "small_strain_j2_elasticity_plane_stress.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainJ2ElasticityPlaneStress::SmallStrainJ2ElasticityPlaneStress()
    : SmallStrainJ2Elasticity3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallStrainJ2ElasticityPlaneStress::SmallStrainJ2ElasticityPlaneStress(const SmallStrainJ2ElasticityPlaneStress& rOther)
    : SmallStrainJ2Elasticity3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainJ2ElasticityPlaneStress::Clone() const
{
    SmallStrainJ2ElasticityPlaneStress::Pointer p_clone(new SmallStrainJ2ElasticityPlaneStress(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainJ2ElasticityPlaneStress::~SmallStrainJ2ElasticityPlaneStress()
{
}

//************************************************************************************
//************************************************************************************

bool& SmallStrainJ2ElasticityPlaneStress::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE)
        rValue = true;

    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void SmallStrainJ2ElasticityPlaneStress::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 4;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ElasticityPlaneStress::CalculateElasticMatrix(VoigtSizeMatrixType& rC, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    // properties should be updated by phase_transformation_system through element before this function is called
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    double NU = rMaterialProperties.GetValue(POISSON_RATIO);
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStress(rC, E, NU);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ElasticityPlaneStress::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
)
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    // properties should be updated by phase_transformation_system through element before this function is called
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    double NU = rMaterialProperties.GetValue(POISSON_RATIO);

    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStress(rStressVector, rStrainVector, E, NU);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ElasticityPlaneStress::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
{
    ConstitutiveLawUtilities<3>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
