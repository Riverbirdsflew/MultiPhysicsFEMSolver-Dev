// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: constitutive_laws_small_strain_application/license.txt
//
//  Main authors:    whf

#pragma once

// System includes

// External includes

// Project includes
#include "small_strain_j2_thermal_plasticity_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/advanced_constitutive_law_small_strain_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SmallStrainJ2ThermalPlasticityPlaneStrain
 * @brief Defines a J2 plasticity constitutive law in 2D (Plane Strain)
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * @note Strain size is 4 (xx, yy, zz, xy).
 * @author whf
 */
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) SmallStrainJ2ThermalPlasticityPlaneStrain
    : public SmallStrainJ2ThermalPlasticity3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo ProcessInfoType;
    typedef SmallStrainJ2ThermalPlasticity3D BaseType;
    typedef std::size_t SizeType;

    using CLutils    = ConstitutiveLawUtilities<3>;
    using AdvCLutils = AdvancedConstitutiveLawSmallStrainUtilities<3>;

    // Counted pointer
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainJ2ThermalPlasticityPlaneStrain);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainJ2ThermalPlasticityPlaneStrain();

    /**
     * @brief Copy constructor.
     */
    SmallStrainJ2ThermalPlasticityPlaneStrain(const SmallStrainJ2ThermalPlasticityPlaneStrain& rOther);

    /**
     * @brief Destructor.
     */
    ~SmallStrainJ2ThermalPlasticityPlaneStrain() override;

    /**
     * @brief Clone function
     * @return A pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return 4;
    };

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

        /**
     * @brief It calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>& rThisVariable,
                           double& rValue) override;

    /**
     * @brief It calculates the value of a specified variable vector
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Vector>& rThisVariable,
                           Vector& rValue) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(
        const ConstitutiveLaw::StrainVectorType &rStrainVector,
        ConstitutiveLaw::StressVectorType &rStressVector,
        ConstitutiveLaw::Parameters &rValues) override;
        

    /**
     * @brief Calculate the thermal strain vector increment
     * @param rThermalStrainVectorIncrement The thermal strain vector increment
     * @param rValues The constitutive law parameters
     * @param rReferenceTemperature The reference temperature
     */
    void CalculateThermalStrainIncrement(Vector& rThermalStrainIncrementVector, ConstitutiveLaw::Parameters &rValues);

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SmallStrainJ2ThermalPlasticityPlaneStrain";
        return buffer.str();
    }
    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Small Strain J2 Thermal Plasticity Plane Strain 2D constitutive law\n";
    };

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the stress and constitutive tensor
     * @param rValues The norm of the deviation stress
     * @param rPlasticStrain
     * @param rAccumulatedPlasticStrain
     */
    void CalculateStressResponse(
            ConstitutiveLaw::Parameters& rValues,
            Vector& rPlasticStrain,
            double& rAccumulatedPlasticStrain,
            Vector& rThermalStrain) override;

    /**
     * @brief This method computes the plastic potential
     * @param DeltaGamma The increment on the Gamma parameter
     * @param NormStressTrial The norm of the stress trial
     * @param YieldFunctionNormalVector The yield function normal vector
     * @param rMaterialProperties The properties of the material
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     */
    void CalculateTangentMatrix(const double DeltaGamma, const double NormStressTrial,
                                const Vector &YieldFunctionNormalVector,
                                ConstitutiveLaw::Parameters& rValues,
                                const double AccumulatedPlasticStrain, Matrix &rElasticityTensor) override;

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor/matrix to be computed
     * @param rMaterialProperties The properties of the material
     */
    void CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues, Matrix &rElasticityTensor) override;

    ///@}
    ///@name Protected  Access
    ///@{
    double& GetThermalStrainNormal() override    { return mThermalStrainNormal; }
    Vector& GetThermalStrain() override          { return mThermalStrain; }
    void SetThermalStrainNormal(const double ToThermalStrainNormal) override {mThermalStrainNormal = ToThermalStrainNormal;}
    void SetThermalStrain(const array_1d<double, 4>& rThermalStrain) {mThermalStrain = rThermalStrain;}
    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class SmallStrainJ2ThermalPlasticityPlaneStrain
} // namespace Kratos.
