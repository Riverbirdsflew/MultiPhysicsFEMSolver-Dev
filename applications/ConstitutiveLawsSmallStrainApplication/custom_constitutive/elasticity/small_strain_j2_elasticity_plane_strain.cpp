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
#include "small_strain_j2_elasticity_plane_strain.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "therm_diff_struct_phase_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

SmallStrainJ2ElasticityPlaneStrain::SmallStrainJ2ElasticityPlaneStrain()
    : SmallStrainJ2Elasticity3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

SmallStrainJ2ElasticityPlaneStrain::SmallStrainJ2ElasticityPlaneStrain(const SmallStrainJ2ElasticityPlaneStrain& rOther)
    : SmallStrainJ2Elasticity3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer SmallStrainJ2ElasticityPlaneStrain::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2ElasticityPlaneStrain>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

SmallStrainJ2ElasticityPlaneStrain::~SmallStrainJ2ElasticityPlaneStrain()
{
}


/***********************************************************************************/
/***********************************************************************************/

Matrix& SmallStrainJ2ElasticityPlaneStrain::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SmallStrainJ2ElasticityPlaneStrain::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& SmallStrainJ2ElasticityPlaneStrain::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void SmallStrainJ2ElasticityPlaneStrain::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW);
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

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ElasticityPlaneStrain::CalculateElasticMatrix(VoigtSizeMatrixType& rC, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    // properties should be updated by phase_transformation_system through element before this function is called
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    double NU = rMaterialProperties.GetValue(POISSON_RATIO);

    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStrain(rC, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ElasticityPlaneStrain::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    // properties should be updated by phase_transformation_system through element before this function is called
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    double NU = rMaterialProperties.GetValue(POISSON_RATIO);

    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2ElasticityPlaneStrain::CalculateCauchyGreenStrain(Parameters& rValues, ConstitutiveLaw::StrainVectorType& rStrainVector)
{
    ConstitutiveLawUtilities<3>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
