// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//            whf

// System includes

// External includes

// Project includes
#include "constitutive_laws_small_strain_application_variables.h"
#include "generic_small_strain_thermal_isotropic_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"

// Yield surfaces
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    if (rElementGeometry.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rElementGeometry.GetValue(REFERENCE_TEMPERATURE);
    } else if (rMaterialProperties.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }
        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();

            // Calculate thermal strain
            Vector ThermalStrainVector = this->GetThermalStrain();
            Vector ThermalStrainIncerementVector(VoigtSize);
            this->CalculateThermalStrainIncrement(ThermalStrainIncerementVector, rValues);
            noalias(ThermalStrainVector) += ThermalStrainIncerementVector;

            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
                //substract themal strain
                if constexpr (Dimension == 2) {
                    noalias(r_strain_vector) -= ThermalStrainVector;
                } else if constexpr (Dimension == 3) {
                    noalias(r_strain_vector) -= ThermalStrainVector;
                }
                noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            } else {
                BaseType::CalculatePK2Stress( r_strain_vector, r_stress_vector, rValues);
                //substract themal strain
                if constexpr (Dimension == 2) {
                    noalias(r_strain_vector) -= ThermalStrainVector;
                } else if constexpr (Dimension == 3) {
                    noalias(r_strain_vector) -= ThermalStrainVector;
                }
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = AdvancedConstitutiveLawSmallStrainUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
        }
        // Calculate thermal strain
        Vector ThermalStrainVector = this->GetThermalStrain();
        Vector ThermalStrainIncerementVector(VoigtSize);
        this->CalculateThermalStrainIncrement(ThermalStrainIncerementVector, rValues);
        noalias(ThermalStrainVector) += ThermalStrainIncerementVector;

        //substract themal strain
        if constexpr (Dimension == 2) {
            noalias(r_strain_vector) -= ThermalStrainVector;
        } else if constexpr (Dimension == 3) {
            noalias(r_strain_vector) -= ThermalStrainVector;
        }

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double threshold = this->GetThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain = this->GetPlasticStrain();

            BoundedArrayType predictive_stress_vector;
            if (r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW)) {
                noalias(predictive_stress_vector) = rValues.GetStressVector();
            } else {
                // S0 = Elastic stress with strain (E-Ep) + S0
                Vector aux_stress = ZeroVector(VoigtSize);
                BaseType::CalculatePK2Stress(r_strain_vector - plastic_strain, aux_stress, rValues);
                this->template AddInitialStressVectorContribution<Vector>(aux_stress);
                noalias(predictive_stress_vector) = aux_stress;
            }

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            array_1d<double, VoigtSize> f_flux; // DF/DS
            array_1d<double, VoigtSize> g_flux; // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment;
            f_flux.clear();
            g_flux.clear();
            plastic_strain_increment.clear();

            // Elastic Matrix
            this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

            // Compute the plastic parameters
            const double F = TConstLawIntegratorType::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain);

            if (F <= std::abs(1.0e-4 * threshold)) { // Elastic case
                noalias(r_integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                TConstLawIntegratorType::IntegrateStressVector(
                    predictive_stress_vector, r_strain_vector, uniaxial_stress,
                    threshold, plastic_denominator, f_flux, g_flux,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, plastic_strain, rValues,
                    characteristic_length);
                noalias(r_integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(rValues, plastic_strain); // this modifies the ConstitutiveMatrix
                }
            }
        }
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    const Vector& rPlasticStrain)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        // Already stored in rValues.GetConstitutiveMatrix()...
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::Secant) {
        const Vector num = prod(rValues.GetConstitutiveMatrix(), rPlasticStrain);
        const double denom = inner_prod(rValues.GetStrainVector(), num);
        noalias(rValues.GetConstitutiveMatrix()) -= outer_prod(num, num) / denom;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 4);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::InitialStiffness) {
        BaseType::CalculateElasticMatrix(rValues.GetConstitutiveMatrix(), rValues);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::OrthogonalSecant) {
        TangentOperatorCalculatorUtility::CalculateOrthogonalSecantTensor(rValues);
    }
}
/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    const double characteristic_length = AdvancedConstitutiveLawSmallStrainUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
    }
    // Calculate thermal strain
    Vector ThermalStrainVector = this->GetThermalStrain();
    Vector ThermalStrainIncerementVector(VoigtSize);
    this->CalculateThermalStrainIncrement(ThermalStrainIncerementVector, rValues);
    noalias(ThermalStrainVector) += ThermalStrainIncerementVector;
    //substract themal strain
    if constexpr (Dimension == 2) {
        noalias(r_strain_vector) -= ThermalStrainVector;
    } else if constexpr (Dimension == 3) {
        noalias(r_strain_vector) -= ThermalStrainVector;
    }

    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // We compute the stress
    // Elastic Matrix
    this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

    // We get some variables
    double threshold = this->GetThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain = this->GetPlasticStrain();

    BoundedArrayType predictive_stress_vector;
    if (r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW)) {
        noalias(predictive_stress_vector) = rValues.GetStressVector();
    } else {
        // Spred = r_constitutive_matrix:(E-Ep) + S0
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);
    }

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    array_1d<double, VoigtSize> f_flux; // DF/DS
    array_1d<double, VoigtSize> g_flux; // DG/DS
    array_1d<double, VoigtSize> plastic_strain_increment;
    f_flux.clear();
    g_flux.clear();
    plastic_strain_increment.clear();

    const double F = TConstLawIntegratorType::CalculatePlasticParameters(
        predictive_stress_vector, r_strain_vector, uniaxial_stress,
        threshold, plastic_denominator, f_flux, g_flux,
        plastic_dissipation, plastic_strain_increment,
        r_constitutive_matrix, rValues, characteristic_length,
        plastic_strain);

    if (F > std::abs(1.0e-4 * threshold)) { // Plastic case
        // While loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, plastic_strain, rValues,
            characteristic_length);
    }

    this->SetPlasticDissipation(plastic_dissipation);
    this->SetPlasticStrain(plastic_strain);
    this->SetThreshold(threshold);
    this->SetThermalStrain(ThermalStrainVector);
    this->SetThermalStrainNormal(ThermalStrainVector[0]);
}


template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::CalculateThermalStrainIncrement(Vector& rThermalStrainIncrementVector, ConstitutiveLaw::Parameters &rValues)
{
    double alpha = CLutils::CalculateMultiPhaseProperties(THERMAL_EXPANSION_COEFFICIENT, rValues);
    const double previous_temperature_gp = AdvCLutils::CalculateInGaussPoint(TEMPERATURE, rValues, 1);
    const double current_temperature_gp = AdvCLutils::CalculateInGaussPoint(TEMPERATURE, rValues, 0);
    const double delta_temperature = current_temperature_gp - previous_temperature_gp;
    if constexpr (Dimension == 2) {
        rThermalStrainIncrementVector[0] = alpha * delta_temperature;
        rThermalStrainIncrementVector[1] = alpha * delta_temperature;
        rThermalStrainIncrementVector[2] = 0.0;
    } else if constexpr (Dimension == 3) {
        rThermalStrainIncrementVector[0] = alpha * delta_temperature;
        rThermalStrainIncrementVector[1] = alpha * delta_temperature;
        rThermalStrainIncrementVector[2] = alpha * delta_temperature;
        rThermalStrainIncrementVector[3] = 0.0;
        rThermalStrainIncrementVector[4] = 0.0;
        rThermalStrainIncrementVector[5] = 0.0;
    }
    // //把一些结果打印到文本文件中
    // const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    // const GeometryType& r_geometry = rValues.GetElementGeometry();
    // const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;
    // if(first_computation){
    //     std::ofstream output_file("thermal_strain.txt", std::ios::trunc);
    //     output_file.close();
    // }
    // if (r_geometry[0].Id() == 1){
    //     const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    //     std::ofstream output_file;
    //     output_file.open("thermal_strain.txt", std::ios::app);
    //     output_file<<"### Previous temperature: " << previous_temperature_gp << " Current temperature: " << current_temperature_gp << " Delta temperature: " << delta_temperature << std::endl;
    //     output_file<<"### Thermal strain increment: " << rThermalStrainIncrementVector << std::endl;
    //     output_file<<"### STEP: " << r_current_process_info[STEP] << std::endl;
    //     output_file<<"### iteration: " << r_current_process_info[NL_ITERATION_NUMBER] << std::endl;
    //     output_file<<"### Thermal expansion coefficient: " << alpha << std::endl;
    //     output_file<<"### Mass Fraction M1: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M1) << std::endl;
    //     output_file<<"### Mass Fraction M2: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M2) << std::endl;
    //     output_file<<"### Mass Fraction M3: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M3) << std::endl;
    //     output_file<<"### Mass Fraction M4: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M4) << std::endl;
    //     output_file<<"### Mass Fraction M5: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M5) << std::endl;
    //     output_file<<"### Mass Fraction M6: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M6) << std::endl;
    //     output_file<<"### Mass Fraction M7: " << r_geometry[0].FastGetSolutionStepValue(MASS_FRACTION_M7) << std::endl;
    //     output_file.close();
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    bool has = false;
    has = BaseType::Has(rThisVariable);
    if (rThisVariable == REFERENCE_TEMPERATURE)
        has = true;
    if (rThisVariable == THERMAL_STRAIN_NORMAL)
        has = true;
    return has;
}

template <class TConstLawIntegratorType>
bool GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == THERMAL_STRAIN_VECTOR) {
        return true;
    }
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (BaseType::Has(rThisVariable))
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    else if (rThisVariable == REFERENCE_TEMPERATURE)
        mReferenceTemperature = rValue;
    else if (rThisVariable == THERMAL_STRAIN_NORMAL)
        mThermalStrainNormal = rValue;
}

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == THERMAL_STRAIN_VECTOR) {
        mThermalStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    rValue = 0.0;
    if (BaseType::Has(rThisVariable))
        BaseType::GetValue(rThisVariable, rValue);
    else if (rThisVariable == REFERENCE_TEMPERATURE)
        rValue = mReferenceTemperature;
    else if (rThisVariable == THERMAL_STRAIN_NORMAL)
        rValue = mThermalStrainNormal;
    return rValue;
}

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == THERMAL_STRAIN_NORMAL) {
        rValue.resize(VoigtSize, false);
        noalias(rValue) = mThermalStrain;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_ERROR_IF_NOT(rElementGeometry[0].SolutionStepsDataHas(TEMPERATURE))  << "The TEMPERATURE variable is not available at the nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT)) << "The THERMAL_EXPANSION_COEFFICIENT is not set in the material properties." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[THERMAL_EXPANSION_COEFFICIENT] < 0.0)   << "The THERMAL_EXPANSION_COEFFICIENT is negative..." << std::endl;
    KRATOS_ERROR_IF_NOT(rElementGeometry.Has(REFERENCE_TEMPERATURE) || rMaterialProperties.Has(REFERENCE_TEMPERATURE)) << "The REFERENCE_TEMPERATURE is not given in the material properties nor via SetValue()" << std::endl;
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;

template class GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
} // namespace Kratos
