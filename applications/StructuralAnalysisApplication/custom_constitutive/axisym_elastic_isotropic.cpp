// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| ANALYSIS
//
//  License:         BSD License
//                   license: StructuralAnalysisApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_utilities/constitutive_law_utilities.h"

#include "structural_analysis_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymElasticIsotropic::AxisymElasticIsotropic()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymElasticIsotropic::AxisymElasticIsotropic(const AxisymElasticIsotropic& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer AxisymElasticIsotropic::Clone() const
{
    AxisymElasticIsotropic::Pointer p_clone(new AxisymElasticIsotropic(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymElasticIsotropic::~AxisymElasticIsotropic()
{
};

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void AxisymElasticIsotropic::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
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

void AxisymElasticIsotropic::CalculateElasticMatrix(VoigtSizeMatrixType& C, ConstitutiveLaw::Parameters& rValues)
{
    //const Properties& MaterialProperties = rValues.GetMaterialProperties();
    //const auto &r_geom = rValues.GetElementGeometry();
    //const auto &r_N = rValues.GetShapeFunctionsValues();
    //const auto &r_process_info = rValues.GetProcessInfo();
    const double E = ConstitutiveLawUtilities<3>::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues);//multiphase, pass through accessor
    const double NU = ConstitutiveLawUtilities<3>::CalculateMultiPhaseProperties(POISSON_RATIO, rValues);//multiphase, pass through accessor
    //const double E = MaterialProperties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);//pass through accessor
    //const double NU = MaterialProperties.GetValue(POISSON_RATIO, r_geom, r_N, r_process_info);//pass through accessor

    if (C.size1() != 4 || C.size2() != 4)
        C.resize(4, 4, false);
    noalias(C) = ZeroMatrix(4, 4);

    const double aux_value = (1.0 - 2.0 * NU);
    const double c1 = E / ((1.0 + NU) * aux_value);
    const double c0 = c1 * (1.0 - NU);
    const double c2 = 0.5 * c1 * aux_value;
    const double c3 = NU * c1;

    C(0, 0) = c0;
    C(0, 1) = c3;
    C(0, 2) = c3;

    C(1, 0) = c3;
    C(1, 1) = c0;
    C(1, 2) = c3;

    C(2, 0) = c3;
    C(2, 1) = c3;
    C(2, 2) = c0;

    C(3, 3) = c2;
}

/***********************************************************************************/
/***********************************************************************************/

void AxisymElasticIsotropic::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    //const Properties& r_material_properties = rValues.GetMaterialProperties();
    //const auto &r_geom = rValues.GetElementGeometry();
    //const auto &r_N = rValues.GetShapeFunctionsValues();
    //const auto &r_process_info = rValues.GetProcessInfo();
    const double E = ConstitutiveLawUtilities<3>::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues);//multiphase, pass through accessor
    const double NU = ConstitutiveLawUtilities<3>::CalculateMultiPhaseProperties(POISSON_RATIO, rValues);//multiphase, pass through accessor
    //const double E  = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);//通过accessor
    //const double NU = r_material_properties.GetValue(POISSON_RATIO, r_geom, r_N, r_process_info);//通过accessor

    const double aux_value = (1.0 - 2.0 * NU);
    const double c1 = E / ((1.0 + NU) * aux_value);
    const double c0 = c1 * (1.0 - NU);
    const double c2 = 0.5 * c1 * aux_value;
    const double c3 = NU * c1;

    rStressVector[0] = (c0 * rStrainVector[0] + c3 * rStrainVector[1] + c3 * rStrainVector[2]);
    rStressVector[1] = (c3 * rStrainVector[0] + c0 * rStrainVector[1] + c3 * rStrainVector[2]);
    rStressVector[2] = (c3 * rStrainVector[0] + c3 * rStrainVector[1] + c0 * rStrainVector[2]);
    rStressVector[3] = (c2 * rStrainVector[3]);
}

//************************************************************************************
//************************************************************************************

void AxisymElasticIsotropic::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveLaw::StrainVectorType& rStrainVector
)
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    const Matrix RightCauchyGreen = prod(trans(F),F);
    rStrainVector[0] = 0.5 * ( RightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( RightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( RightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = RightCauchyGreen( 0, 1 );
}

} // Namespace Kratos
