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

# pragma once

// System includes

// External includes

// Project includes
#include "small_strain_j2_thermal_elasticity_3d.h"

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
 * @class SmallStrainJ2ThermalElasticityPlaneStrain
 * @ingroup constitutive_laws_small_strain_application
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 */
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) SmallStrainJ2ThermalElasticityPlaneStrain
    : public SmallStrainJ2Elasticity3D
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = SmallStrainJ2Elasticity3D;

        /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;

    /// Counted pointer of SmallStrainJ2ThermalElasticityPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainJ2ThermalElasticityPlaneStrain);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainJ2ThermalElasticityPlaneStrain() 
    {}

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SmallStrainJ2ThermalElasticityPlaneStrain>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~SmallStrainJ2ThermalElasticityPlaneStrain() override {}

    /**
    * Copy constructor.
    */
    SmallStrainJ2ThermalElasticityPlaneStrain(const SmallStrainJ2ThermalElasticityPlaneStrain &rOther)
        : BaseType(rOther)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() const override
    {
        return VoigtSize;
    }

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
    * @brief It calculates the constitutive matrix rConstitutiveMatrix
    * @param rConstitutiveMatrix The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues
        ) override;

    /**
     * @brief Calculate the thermal strain vector increment
     * @param rThermalStrainVectorIncrement The thermal strain vector increment
     * @param rValues The constitutive law parameters
     * @param rReferenceTemperature The reference temperature
     */
    virtual void CalculateThermalStrainIncrement(Vector& rThermalStrainIncrementVector, ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief It calculates the strain vector
     * @param rValues The internal values of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw::StrainVectorType& rStrainVector
        ) override;

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

    ///@}
    ///@name Protected Operations
    ///@{

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainJ2Elasticity3D)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainJ2Elasticity3D)
    }


}; // Class SmallStrainJ2ThermalElasticityPlaneStrain
}  // namespace Kratos.
