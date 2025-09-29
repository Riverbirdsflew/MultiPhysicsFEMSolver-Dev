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
#include "small_strain_j2_thermal_plasticity_3d.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "structural_analysis_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainJ2ThermalPlasticity3D::SmallStrainJ2ThermalPlasticity3D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainJ2ThermalPlasticity3D::SmallStrainJ2ThermalPlasticity3D(const SmallStrainJ2ThermalPlasticity3D &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainJ2ThermalPlasticity3D::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2ThermalPlasticity3D>(SmallStrainJ2ThermalPlasticity3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainJ2ThermalPlasticity3D::~SmallStrainJ2ThermalPlasticity3D()
{
}

//************************************************************************************
//************************************************************************************

bool SmallStrainJ2ThermalPlasticity3D::Has(const Variable<double>& rThisVariable)
{
    bool has = false;
    has = BaseType::Has(rThisVariable);
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        has = true;
    }
    if(rThisVariable == THERMAL_STRAIN_NORMAL){
        has = true;
    }
    return has;
}

bool SmallStrainJ2ThermalPlasticity3D::Has(const Variable<Vector>& rThisVariable)
{
    bool has = false;
    has = BaseType::Has(rThisVariable);
    if(rThisVariable == THERMAL_STRAIN_VECTOR){
        has = true;
    }
    return has;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (BaseType::Has(rThisVariable)){
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
    else if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        mAccumulatedPlasticStrain = rValue;
    }
    else if(rThisVariable == THERMAL_STRAIN_NORMAL){
        mThermalStrainNormal = rValue;
    }
}

void SmallStrainJ2ThermalPlasticity3D::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    
    if (BaseType::Has(rThisVariable)){
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
    else if(rThisVariable == THERMAL_STRAIN_VECTOR){
        mThermalStrain = rValue;
    }
}

//************************************************************************************
//************************************************************************************

double& SmallStrainJ2ThermalPlasticity3D::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (BaseType::Has(rThisVariable)){
        BaseType::GetValue(rThisVariable, rValue);
    }
    else if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }
    else if(rThisVariable == THERMAL_STRAIN_NORMAL){
        rValue = mThermalStrainNormal;
    }
    return rValue;
}

Vector& SmallStrainJ2ThermalPlasticity3D::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (BaseType::Has(rThisVariable)){
        BaseType::GetValue(rThisVariable, rValue);
    }
    else if(rThisVariable == THERMAL_STRAIN_VECTOR){
        rValue.resize(6, false);
        noalias(rValue) = mThermalStrain;
    }
    return rValue;
}
//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mPlasticStrain = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrain = 0.0;
    mThermalStrain = ZeroVector(this->GetStrainSize());
    mThermalStrainNormal = 0.0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::InitializeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::InitializeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::InitializeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::InitializeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::FinalizeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    Vector plastic_strain;
    double accumulated_plastic_strain;
    Vector thermal_strain;
    this->CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain, thermal_strain);
    mPlasticStrain = plastic_strain;
    mAccumulatedPlasticStrain = accumulated_plastic_strain;
    this->SetThermalStrain(thermal_strain);
    this->SetThermalStrainNormal(thermal_strain[0]);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::FinalizeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::FinalizeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::FinalizeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // In small deformation is the same as compute Cauchy
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // In small deformation is the same as compute Cauchy
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    Vector plastic_strain;
    Vector thermal_strain;
    double accumulated_plastic_strain;

    this->CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain, thermal_strain);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rPlasticStrain,
    double& rAccumulatedPlasticStrain,
    Vector& rThermalStrain
    )
{
    rPlasticStrain.resize(6, false);
    rThermalStrain.resize(6, false);
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
    Vector ThermalStrainIncerementVector(6);
    this->CalculateThermalStrainIncrement(ThermalStrainIncerementVector, rValues);
    noalias(rThermalStrain) += ThermalStrainIncerementVector;

    //substract themal strain
    noalias(r_strain_vector) -= rThermalStrain;

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        const double mu = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (2. * (1. + CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
        const double bulk_modulus = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (3. * (1. - 2. * CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
        double trial_yield_function;

        rPlasticStrain = mPlasticStrain;
        rAccumulatedPlasticStrain = mAccumulatedPlasticStrain;

        Matrix elastic_tensor;
        elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(rValues, elastic_tensor);
        Vector trial_stress(6);
        noalias(trial_stress) = prod(elastic_tensor, r_strain_vector - mPlasticStrain);

        // trial_stress_dev = sigma - 1/3 tr(sigma) * I
        Vector trial_stress_dev = trial_stress;
        const double trace = 1. / 3. * (trial_stress(0) + trial_stress(1) + trial_stress(2));
        trial_stress_dev(0) -= trace;
        trial_stress_dev(1) -= trace;
        trial_stress_dev(2) -= trace;
        const double norm_trial_stress_dev = std::sqrt(trial_stress_dev(0) * trial_stress_dev(0) +
                                        trial_stress_dev(1) * trial_stress_dev(1) +
                                        trial_stress_dev(2) * trial_stress_dev(2) +
                                        2. * trial_stress_dev(3) * trial_stress_dev(3) +
                                        2. * trial_stress_dev(4) * trial_stress_dev(4) +
                                        2. * trial_stress_dev(5) * trial_stress_dev(5));
                                        
        trial_yield_function = this->YieldFunction(norm_trial_stress_dev, rValues, mAccumulatedPlasticStrain, 0.0);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector = trial_stress;
            }
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                r_constitutive_matrix = elastic_tensor;
            }
        } else {
            // INELASTIC
            double D_gamma = 0;
            Vector yield_function_normal_vector = trial_stress_dev / norm_trial_stress_dev;
            // Get the Dgammma
            D_gamma = GetDgamma(norm_trial_stress_dev, rValues, mAccumulatedPlasticStrain);
            
            // We update the stress
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector(0) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    trial_stress_dev(0) - 2. * mu * D_gamma * yield_function_normal_vector(0);
                r_stress_vector(1) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    trial_stress_dev(1) - 2. * mu * D_gamma * yield_function_normal_vector(1);
                r_stress_vector(2) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    trial_stress_dev(2) - 2. * mu * D_gamma * yield_function_normal_vector(2);
                r_stress_vector(3) =
                    trial_stress_dev(3) - 2. * mu * D_gamma * yield_function_normal_vector(3);
                r_stress_vector(4) =
                    trial_stress_dev(4) - 2. * mu * D_gamma * yield_function_normal_vector(4);
                r_stress_vector(5) =
                    trial_stress_dev(5) - 2. * mu * D_gamma * yield_function_normal_vector(5);
            }

            rPlasticStrain(0) += D_gamma * yield_function_normal_vector(0);
            rPlasticStrain(1) += D_gamma * yield_function_normal_vector(1);
            rPlasticStrain(2) += D_gamma * yield_function_normal_vector(2);
            rPlasticStrain(3) += D_gamma * yield_function_normal_vector(3) * 2;
            rPlasticStrain(4) += D_gamma * yield_function_normal_vector(4) * 2;
            rPlasticStrain(5) += D_gamma * yield_function_normal_vector(5) * 2;
            rAccumulatedPlasticStrain += Globals::SqrtTwoThird * D_gamma;

            // We update the tangent tensor
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                CalculateTangentMatrix(D_gamma, norm_trial_stress_dev, yield_function_normal_vector,
                                       rValues, rAccumulatedPlasticStrain, r_constitutive_matrix);
            }
        }
    }
}

//************************************************************************************
//************************************************************************************

double& SmallStrainJ2ThermalPlasticity3D::CalculateValue(
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
                                        2. * stress_dev(3) * stress_dev(3) +
                                        2. * stress_dev(4) * stress_dev(4) +
                                        2. * stress_dev(5) * stress_dev(5));
    }   
    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainJ2ThermalPlasticity3D::CalculateValue(
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
        
        if(rValue.size() != 6){
            rValue.resize(6, false);
        }
        this->CalculateMaterialResponseCauchy(rParameterValues);
        noalias(rValue) = rParameterValues.GetStressVector();
        // Previous flags restored
    }
    return(rValue);
}


//************************************************************************************
//************************************************************************************

double SmallStrainJ2ThermalPlasticity3D::GetYieldStress(ConstitutiveLaw::Parameters& rValues,
    const double accumulated_plastic_strain)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    KRATOS_ERROR_IF(!(r_material_properties.HasM4dTable(2,TEMPERATURE,
                                                        PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT,
                                                        PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT,
                                                        PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT, 
                                                        EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT))) 
                    << "No M4dTable for EqStress - EqPlasticStrain and temperature " << std::endl;

    const M4dTable plastic_strain_stress_table = r_material_properties.GetM4dTable(2,TEMPERATURE,
                                                                             PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT,
                                                                             PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT,
                                                                             PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT, 
                                                                             EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT);

    const double eq_plastic_strain = accumulated_plastic_strain;
    const double current_temperature = AdvancedConstitutiveLawSmallStrainUtilities<6>::CalculateInGaussPoint(TEMPERATURE, rValues, 0);

    const double yield_stress = plastic_strain_stress_table.GetInterpolateValue({current_temperature,eq_plastic_strain,eq_plastic_strain,eq_plastic_strain});
    return yield_stress;
}


//************************************************************************************
//************************************************************************************

double SmallStrainJ2ThermalPlasticity3D::GetDgamma(
    const double NormTrialStressDev,
    ConstitutiveLaw::Parameters& rValues,
    const double AccumulatedPlasticStrainOld
)
{
    const double mu = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (2. * (1. + CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
    double yf_new = this->GetYieldStress(rValues, AccumulatedPlasticStrainOld);
    const double tolerance = 1e-13 * yf_new;
    double accumulated_plastic_strain = AccumulatedPlasticStrainOld;
    double Dgamma = 0.0;
    SizeType num_iterations = 0;

    while (yf_new > tolerance && num_iterations < 300)
    {
        const double yf = this->YieldFunction(NormTrialStressDev, rValues, accumulated_plastic_strain, Dgamma);
        const double dyf = -2. * mu * (1. + this->GetHardeningModulus(rValues, accumulated_plastic_strain) / (3. * mu));
        Dgamma -= yf / dyf;
        accumulated_plastic_strain = AccumulatedPlasticStrainOld + Globals::SqrtTwoThird * Dgamma;
        yf_new = std::abs(yf);
        num_iterations++;
    }
    // limit the number of iterations and throw a warning
    if (num_iterations >= 300)
    {
        KRATOS_WARNING(this->Info()) << "No convergence in the plastic strain rate calculation" << std::endl;
    }

    return Dgamma;
}

//************************************************************************************
//************************************************************************************

double SmallStrainJ2ThermalPlasticity3D::YieldFunction(
    const double NormTrialStressDev,
    ConstitutiveLaw::Parameters& rValues,
    const double AccumulatedPlasticStrain,
    const double Dgamma
    )
{
    const double mu = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (2. * (1. + CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));

    return (NormTrialStressDev - Globals::SqrtTwoThird * this->GetYieldStress(rValues, AccumulatedPlasticStrain) - 2. * mu * Dgamma);
}

//************************************************************************************
//************************************************************************************
double SmallStrainJ2ThermalPlasticity3D::GetHardeningModulus(
    ConstitutiveLaw::Parameters& rValues,
    const double AccumulatedPlasticStrain
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto& plastic_strain_stress_table = r_material_properties.GetM4dTable(2,TEMPERATURE,
                                                                PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT,
                                                                PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT,
                                                                PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT, 
                                                                EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT);
    
    KRATOS_ERROR_IF(!(plastic_strain_stress_table.Size(0)) && (plastic_strain_stress_table.Size(1))) 
                        << "No M4dTable for EqStress - EqPlasticStrain and temperature " << std::endl;
    
    double current_temperature = AdvancedConstitutiveLawSmallStrainUtilities<6>::CalculateInGaussPoint(TEMPERATURE, rValues, 0);
    double eq_plastic_strain = AccumulatedPlasticStrain;
    double hardening_modulus = plastic_strain_stress_table.GetDerivativeValue(1,{current_temperature,eq_plastic_strain,eq_plastic_strain,eq_plastic_strain});
    //std::cout << "##### hardening_modulus: " << hardening_modulus << std::endl;
    return hardening_modulus;
}


//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateElasticMatrix(
    ConstitutiveLaw::Parameters& rValues, Matrix &rElasticMatrix)
{
    const double E = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues);
    const double poisson_ratio = CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues);
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. * (1. + poisson_ratio));

    if (rElasticMatrix.size1() != 6 || rElasticMatrix.size2() != 6)
        rElasticMatrix.resize(6, 6, false);
    rElasticMatrix.clear();

    rElasticMatrix(0, 0) = lambda + 2. * mu;
    rElasticMatrix(0, 1) = lambda;
    rElasticMatrix(0, 2) = lambda;
    rElasticMatrix(1, 0) = lambda;
    rElasticMatrix(1, 1) = lambda + 2. * mu;
    rElasticMatrix(1, 2) = lambda;
    rElasticMatrix(2, 0) = lambda;
    rElasticMatrix(2, 1) = lambda;
    rElasticMatrix(2, 2) = lambda + 2. * mu;
    rElasticMatrix(3, 3) = mu;
    rElasticMatrix(4, 4) = mu;
    rElasticMatrix(5, 5) = mu;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::CalculateTangentMatrix(
        const double Dgamma, const double NormTrialStressDev,
        const Vector &rYFNormalVector,
        ConstitutiveLaw::Parameters &rValues,
        const double AccumulatedPlasticStrain,
        Matrix &rTMatrix)
{
    const double hardening_modulus = this->GetHardeningModulus(rValues, AccumulatedPlasticStrain);
    const double mu = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (2. * (1. + CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));
    const double bulk_modulus = CLutils::CalculateMultiPhaseProperties(YOUNG_MODULUS, rValues) / (3. * (1. - 2. * CLutils::CalculateMultiPhaseProperties(POISSON_RATIO, rValues)));

    const double theta_new = 1 - (2. * mu * Dgamma) / NormTrialStressDev;
    const double theta_new_b = 1. / (1. + hardening_modulus / (3. * mu)) - (1. - theta_new);

    rTMatrix(0, 0) = bulk_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(0)));
    rTMatrix(0, 1) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(1)));
    rTMatrix(0, 2) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(2)));
    rTMatrix(0, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(3)));
    rTMatrix(0, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(4)));
    rTMatrix(0, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(5)));

    rTMatrix(1, 0) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(0)));
    rTMatrix(1, 1) = bulk_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(1)));
    rTMatrix(1, 2) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(2)));
    rTMatrix(1, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(3)));
    rTMatrix(1, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(4)));
    rTMatrix(1, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(5)));

    rTMatrix(2, 0) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(0)));
    rTMatrix(2, 1) = bulk_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(1)));
    rTMatrix(2, 2) = bulk_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(2)));
    rTMatrix(2, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(3)));
    rTMatrix(2, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(4)));
    rTMatrix(2, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(5)));

    rTMatrix(3, 0) = -(2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(0)));
    rTMatrix(3, 1) = -(2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(1)));
    rTMatrix(3, 2) = -(2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(2)));
    rTMatrix(3, 3) = mu * theta_new - (2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(3)));
    rTMatrix(3, 4) = -(2. *mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(4)));
    rTMatrix(3, 5) = -(2. *mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(5)));

    rTMatrix(4, 0) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(0)));
    rTMatrix(4, 1) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(1)));
    rTMatrix(4, 2) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(2)));
    rTMatrix(4, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(3)));
    rTMatrix(4, 4) = mu * theta_new - (2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(4)));
    rTMatrix(4, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(5)));

    rTMatrix(5, 0) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(0)));
    rTMatrix(5, 1) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(1)));
    rTMatrix(5, 2) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(2)));
    rTMatrix(5, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(3)));
    rTMatrix(5, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(4)));
    rTMatrix(5, 5) = mu * theta_new - (2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(5)));
}

//************************************************************************************
//************************************************************************************
void SmallStrainJ2ThermalPlasticity3D::CalculateThermalStrainIncrement(
    Vector& rThermalStrainIncrementVector, 
    ConstitutiveLaw::Parameters &rValues
    )
{
    const double previous_temperature_gp = AdvCLutils::CalculateInGaussPoint(TEMPERATURE, rValues, 1);
    const double current_temperature_gp = AdvCLutils::CalculateInGaussPoint(TEMPERATURE, rValues, 0);
    const double delta_temperature = current_temperature_gp - previous_temperature_gp;
    if (this->WorkingSpaceDimension() == 3) {
        rThermalStrainIncrementVector[0] = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues) * delta_temperature;
        rThermalStrainIncrementVector[1] = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues) * delta_temperature;
        rThermalStrainIncrementVector[2] = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues) * delta_temperature;
        rThermalStrainIncrementVector[3] = 0.0;
        rThermalStrainIncrementVector[4] = 0.0;
        rThermalStrainIncrementVector[5] = 0.0;
    } else if (this->WorkingSpaceDimension() == 2) {
        rThermalStrainIncrementVector[0] = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues) * delta_temperature;
        rThermalStrainIncrementVector[1] = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues) * delta_temperature;
        rThermalStrainIncrementVector[2] = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues) * delta_temperature;;
        rThermalStrainIncrementVector[3] = 0.0;
    }
    else{
        KRATOS_ERROR << "The dimension of the problem is not correct" << std::endl;
    }
}
//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->GetStrainSize();
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************

int SmallStrainJ2ThermalPlasticity3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT));

    return 0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.save("mThermalStrainNormal", mThermalStrainNormal);
    rSerializer.save("mThermalStrain", mThermalStrain);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2ThermalPlasticity3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.load("mThermalStrainNormal", mThermalStrainNormal);
    rSerializer.load("mThermalStrain", mThermalStrain);
}

} /* namespace Kratos.*/
