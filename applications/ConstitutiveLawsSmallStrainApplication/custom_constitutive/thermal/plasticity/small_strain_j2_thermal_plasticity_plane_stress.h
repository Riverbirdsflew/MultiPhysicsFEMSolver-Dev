// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: constitutive_laws_small_strain_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

# pragma once

// System includes

// External includes

// Project includes
#include "small_strain_j2_thermal_plasticity_3d.h"

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
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 */
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) SmallStrainJ2ThermalPlasticityPlaneStrain
    : public SmallStrainJ2ThermalPlasticity3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef SmallStrainJ2ThermalPlasticity3D BaseType;

    /// Counted pointer of LinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainJ2ThermalPlasticityPlaneStrain);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainJ2ThermalPlasticityPlaneStrain() 
    {}

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SmallStrainJ2ThermalPlasticityPlaneStrain>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~SmallStrainJ2ThermalPlasticityPlaneStrain() override {}

    /**
    * Copy constructor.
    */
    SmallStrainJ2ThermalPlasticityPlaneStrain(const SmallStrainJ2ThermalPlasticityPlaneStrain &rOther)
        : BaseType(rOther)
    {
    }

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

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
    ///@}
    ///@name Access
    ///@{


    /**
     * @brief Calculate the thermal strain vector increment
     * @param rThermalStrainVectorIncrement The thermal strain vector increment
     * @param rValues The constitutive law parameters
     * @param rReferenceTemperature The reference temperature
     */
    void CalculateThermalStrainIncrement(Vector& rThermalStrainIncrementVector, ConstitutiveLaw::Parameters &rValues);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
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

    double& GetThermalStrainNormal() override    { return mThermalStrainNormal; }
    Vector& GetThermalStrain() override          { return mThermalStrain; }
    void SetThermalStrainNormal(const double ToThermalStrainNormal) override {mThermalStrainNormal = ToThermalStrainNormal;}
    void SetThermalStrain(const array_1d<double, 4>& rThermalStrain) {mThermalStrain = rThermalStrain;}

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainJ2ThermalPlasticity3D)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainJ2ThermalPlasticity3D)
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
