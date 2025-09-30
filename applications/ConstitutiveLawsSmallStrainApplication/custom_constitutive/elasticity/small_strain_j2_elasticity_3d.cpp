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
#include "small_strain_j2_elasticity_3d.h"
#include "includes/checks.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

SmallStrainJ2Elasticity3D::SmallStrainJ2Elasticity3D()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

SmallStrainJ2Elasticity3D::SmallStrainJ2Elasticity3D(const SmallStrainJ2Elasticity3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer SmallStrainJ2Elasticity3D::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2Elasticity3D>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

SmallStrainJ2Elasticity3D::~SmallStrainJ2Elasticity3D()
{
}

/***********************************************************************************/
/***********************************************************************************/


// NOte: for small strains case, CalculateMaterialResponseCauchy should be called 

void  SmallStrainJ2Elasticity3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void SmallStrainJ2Elasticity3D::CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    Flags &r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType &r_strain_vector = rValues.GetStrainVector();

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since in small strains, any strain measure works, we can use CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
        AddInitialStressVectorContribution<StressVectorType>(r_stress_vector);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

double& SmallStrainJ2Elasticity3D::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    Flags &r_constitutive_law_options = rParameterValues.GetOptions();
    ConstitutiveLaw::StrainVectorType &r_strain_vector = rParameterValues.GetStrainVector();
    ConstitutiveLaw::StressVectorType &r_stress_vector = rParameterValues.GetStressVector();

    if (rThisVariable == STRAIN_ENERGY) {
        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        }
        this->CalculatePK2Stress(r_strain_vector, r_stress_vector, rParameterValues);
        rValue = 0.5 * inner_prod(r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SmallStrainJ2Elasticity3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN) {

        Flags& r_flags = rParameterValues.GetOptions();

        ConstitutiveLaw::StrainVectorType& r_strain_vector = rParameterValues.GetStrainVector();

        if( r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
            //Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
            CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        }
        AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = r_strain_vector;

    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR) {

        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        // Set flags to only compute the stress
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        this->CalculateMaterialResponsePK2(rParameterValues);
        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
        return rValue;

    } else if (rThisVariable == INITIAL_STRAIN_VECTOR) {
        if (this->HasInitialState()) {
	    if (rValue.size() != GetStrainSize()) {
	        rValue.resize(GetStrainSize());
	    }
	    noalias(rValue) = GetInitialState().GetInitialStrainVector();
        } else {
            noalias(rValue) = ZeroVector(0);
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& SmallStrainJ2Elasticity3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CONSTITUTIVE_MATRIX) {
        this->CalculateElasticMatrix(rValue, rParameterValues);
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 6;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/

int SmallStrainJ2Elasticity3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] < 0.0) << "YOUNG_MODULUS is negative." << std::endl;

    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

std::string SmallStrainJ2Elasticity3D::Info() const
{
    return "SmallStrainJ2Elasticity3D ConstitutiveLaw instance";
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::CheckClearElasticMatrix(ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix)
{
    const SizeType size_system = this->GetStrainSize();
    if (rConstitutiveMatrix.size1() != size_system || rConstitutiveMatrix.size2() != size_system)
        rConstitutiveMatrix.resize(size_system, size_system, false);
    rConstitutiveMatrix.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    // properties should be updated by phase_transformation_system through element before this function is called
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    double NU = rMaterialProperties.GetValue(POISSON_RATIO);
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(rConstitutiveMatrix, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    // properties should be updated by phase_transformation_system through element before this function is called
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS);
    double NU = rMaterialProperties.GetValue(POISSON_RATIO);

    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallStrainJ2Elasticity3D::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveLaw::StrainVectorType& rStrainVector
    )
{
    ConstitutiveLawUtilities<6>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
