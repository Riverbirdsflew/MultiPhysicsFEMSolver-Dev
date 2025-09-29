// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Alfredo Huespe
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "small_strain_j2_thermal_plasticity_plane_strain_2d.h"
#include "constitutive_laws_small_strain_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainJ2ThermalPlasticityPlaneStrain2D::SmallStrainJ2ThermalPlasticityPlaneStrain2D()
    : SmallStrainJ2ThermalPlasticity3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainJ2ThermalPlasticityPlaneStrain2D::SmallStrainJ2ThermalPlasticityPlaneStrain2D(const SmallStrainJ2ThermalPlasticityPlaneStrain2D &rOther)
    : SmallStrainJ2ThermalPlasticity3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainJ2ThermalPlasticityPlaneStrain2D::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2ThermalPlasticityPlaneStrain2D>(SmallStrainJ2ThermalPlasticityPlaneStrain2D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainJ2ThermalPlasticityPlaneStrain2D::~SmallStrainJ2ThermalPlasticityPlaneStrain2D()
{
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainJ2ThermalPlasticityPlaneStrain2D::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (BaseType::Has(rThisVariable)){
        BaseType::GetValue(rThisVariable, rValue);
    }
    else if(rThisVariable == THERMAL_STRAIN_VECTOR){
        rValue.resize(4, false);
        noalias(rValue) = mThermalStrain;
    }
    return rValue;
}
//************************************************************************************
//************************************************************************************

double& SmallStrainJ2ThermalPlasticityPlaneStrain2D::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN){

        rValue = Globals::SqrtTwoThird * std::sqrt(inner_prod(mPlasticStrain, mPlasticStrain));
        //rValue = this->GetValue(ACCUMULATED_PLASTIC_STRAIN, rValue);
    }else if(rThisVariable == UNIAXIAL_STRESS){
        this->CalculateMaterialResponseCauchy(rValues);
        const Vector& r_stress_vector = rValues.GetStressVector();
        Vector stress_dev = r_stress_vector;
        const double I1 = 1. / 3. * (r_stress_vector(0) + r_stress_vector(1) + r_stress_vector(2));
        stress_dev(0) -= I1;
        stress_dev(1) -= I1;
        stress_dev(2) -= I1;

        rValue = std::sqrt(3.0/2.0) * std::sqrt(stress_dev(0) * stress_dev(0) +
                                        stress_dev(1) * stress_dev(1) +
                                        stress_dev(2) * stress_dev(2) +
                                        2. * stress_dev(3) * stress_dev(3));
    }   
    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainJ2ThermalPlasticityPlaneStrain2D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN) {
        const SizeType space_dimension = this->WorkingSpaceDimension();

        // Compute total deformation gradient
        const Matrix& r_F = rParameterValues.GetDeformationGradientF();
        KRATOS_DEBUG_ERROR_IF(r_F.size1()!= space_dimension || r_F.size2() != space_dimension)
            << "expected size of F " << space_dimension << "x" << space_dimension
            << ", got " << r_F.size1() << "x" << r_F.size2() << std::endl;

        const Matrix C_tensor = prod(trans(r_F),r_F);
        CLutils::CalculateGreenLagrangianStrain(C_tensor, rValue);
    }
    if(rThisVariable == STRESSES){
        
        if(rValue.size() != 4){
            rValue.resize(4, false);
        }
        this->CalculateMaterialResponseCauchy(rParameterValues);
        noalias(rValue) = rParameterValues.GetStressVector();
        // Previous flags restored
    }
    return(rValue);
}
//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticityPlaneStrain2D::CalculateStressResponse(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rPlasticStrain,
        double& rAccumulatedPlasticStrain,
        Vector& rThermalStrain)
{
    rPlasticStrain.resize(4, false);
    rThermalStrain.resize(4, false);
    Flags& r_constitutive_law_options = rValues.GetOptions();
    Vector& r_strain_vector = rValues.GetStrainVector();
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    if( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        //this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Calculate thermal strain
    rThermalStrain = this->GetThermalStrain();
    Vector ThermalStrainIncerementVector(4);
    this->CalculateThermalStrainIncrement(ThermalStrainIncerementVector, rValues);
    noalias(rThermalStrain) += ThermalStrainIncerementVector;

    //substract themal strain
    noalias(r_strain_vector) -= rThermalStrain;

    if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();

        const double mu = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (2. * (1. + CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
        const double bulk_modulus = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (3. * (1. - 2. * CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
        const double sqrt_two_thirds = std::sqrt(2. / 3.); // = 0.8164965809277260
        double trial_yield_function;

        rPlasticStrain = mPlasticStrain;
        rAccumulatedPlasticStrain = mAccumulatedPlasticStrain;

        Matrix elastic_tensor;
        elastic_tensor.resize(4, 4, false);
        CalculateElasticMatrix(rValues, elastic_tensor);
        Vector trial_stress(4);
        noalias(trial_stress) = prod(elastic_tensor, r_strain_vector - mPlasticStrain);

        // stress_trial_dev = sigma - 1/3 tr(sigma) * I
        Vector stress_trial_dev = trial_stress;

        const double trace = 1. / 3. * (trial_stress(0) + trial_stress(1) + trial_stress(2));
        stress_trial_dev(0) -= trace;
        stress_trial_dev(1) -= trace;
        stress_trial_dev(2) -= trace;
        const double norm_dev_stress = std::sqrt(stress_trial_dev(0) * stress_trial_dev(0) +
                                           stress_trial_dev(1) * stress_trial_dev(1) +
                                           stress_trial_dev(2) * stress_trial_dev(2) +
                                           2. * stress_trial_dev(3) * stress_trial_dev(3));
        trial_yield_function = this->YieldFunction(norm_dev_stress, rValues, mAccumulatedPlasticStrain, 0.0);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector = trial_stress;
            }
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                tangent_tensor = elastic_tensor;
            }
        } else {
            // INELASTIC
            double D_gamma = 0;
            Vector yield_function_normal_vector = stress_trial_dev / norm_dev_stress;

            D_gamma = GetDgamma(norm_dev_stress, rValues, mAccumulatedPlasticStrain);

            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector(0) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    stress_trial_dev(0) - 2. * mu * D_gamma * yield_function_normal_vector(0);
                r_stress_vector(1) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    stress_trial_dev(1) - 2. * mu * D_gamma * yield_function_normal_vector(1);
                r_stress_vector(2) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    stress_trial_dev(2) - 2. * mu * D_gamma * yield_function_normal_vector(2);
                r_stress_vector(3) =
                        stress_trial_dev(3) - 2. * mu * D_gamma * yield_function_normal_vector(3);
                }

            rPlasticStrain(0) += D_gamma * yield_function_normal_vector(0);
            rPlasticStrain(1) += D_gamma * yield_function_normal_vector(1);
            rPlasticStrain(2) += D_gamma * yield_function_normal_vector(2);
            rPlasticStrain(3) += D_gamma * yield_function_normal_vector(3) * 2;
            rAccumulatedPlasticStrain += sqrt_two_thirds * D_gamma;

            // We update the tangent tensor
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                CalculateTangentMatrix(D_gamma, norm_dev_stress, yield_function_normal_vector,
                                       rValues, rAccumulatedPlasticStrain, tangent_tensor);
            }
        }
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticityPlaneStrain2D::CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues, Matrix &rElasticityTensor)
{
    const double E = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues);
    const double poisson_ratio = CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues);
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticityTensor.size1() != 4 || rElasticityTensor.size2() != 4)
        rElasticityTensor.resize(4, 4, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2. * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(0, 3) = 0.;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2. * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(1, 3) = 0.;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2. * mu;
    rElasticityTensor(2, 3) = 0.;
    rElasticityTensor(3, 0) = 0.;
    rElasticityTensor(3, 1) = 0.;
    rElasticityTensor(3, 2) = 0.;
    rElasticityTensor(3, 3) = mu;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticityPlaneStrain2D::CalculateTangentMatrix(const double DeltaGamma, const double NormStressTrial,
                                                             const Vector &YieldFunctionNormalVector,
                                                             ConstitutiveLaw::Parameters &rValues,
                                                             const double AccumulatedPlasticStrain,
                                                             Matrix &rElasticityTensor)
{
    const double hardening_modulus = this->GetHardeningModulus(rValues, AccumulatedPlasticStrain);
    const double mu = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (2. * (1. + CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
    const double bulk_modulus = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (3. * (1. - 2. * CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));

    const double theta_new = 1 - (2. * mu * DeltaGamma) / NormStressTrial;
    const double theta_new_b = 1. / (1. + hardening_modulus / (3. * mu)) - (1. - theta_new);

    rElasticityTensor(0, 0) = bulk_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(0)));
    rElasticityTensor(0, 1) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(1)));
    rElasticityTensor(0, 2) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(2)));
    rElasticityTensor(0, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(3)));

    rElasticityTensor(1, 0) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(0)));
    rElasticityTensor(1, 1) = bulk_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(1)));
    rElasticityTensor(1, 2) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(2)));
    rElasticityTensor(1, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(3)));

    rElasticityTensor(2, 0) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(0)));
    rElasticityTensor(2, 1) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(1)));
    rElasticityTensor(2, 2) = bulk_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(2)));
    rElasticityTensor(2, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(3)));

    rElasticityTensor(3, 0) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(0)));
    rElasticityTensor(3, 1) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(1)));
    rElasticityTensor(3, 2) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(2)));
    rElasticityTensor(3, 3) = mu * theta_new - (2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(3)));
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticityPlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 4;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticityPlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainJ2ThermalPlasticity3D);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticityPlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainJ2ThermalPlasticity3D);
}

} /* namespace Kratos.*/
