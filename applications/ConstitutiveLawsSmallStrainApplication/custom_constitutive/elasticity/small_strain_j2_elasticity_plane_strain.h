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

#pragma once

// System includes

// External includes

// Project includes
#include "small_strain_j2_elasticity_3d.h"

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
 * @class SmallStrainJ2ElasticityPlaneStrain
 * @ingroup constitutive_laws_small_strain_application
 * @brief This class defines a small deformation linear elastic constitutive model for plane strain cases
 * @details This class derives from the linear elastic case on 3D
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) SmallStrainJ2ElasticityPlaneStrain
    : public SmallStrainJ2Elasticity3D
{
public:
    ///@name Type Definitions
    ///@{

    /// The process info definition
    typedef ProcessInfo      ProcessInfoType;

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw       CLBaseType;

    /// The base class SmallStrainJ2Elasticity3D type definition
    typedef SmallStrainJ2Elasticity3D      BaseType;

    /// The size type definition
    typedef std::size_t             SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 4;

    /// Counted pointer of SmallStrainJ2ElasticityPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION( SmallStrainJ2ElasticityPlaneStrain );

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainJ2ElasticityPlaneStrain();

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    SmallStrainJ2ElasticityPlaneStrain (const SmallStrainJ2ElasticityPlaneStrain& rOther);


    /**
     * @brief Destructor.
     */
    ~SmallStrainJ2ElasticityPlaneStrain() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

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

    /**
     * @brief  Itreturns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */

    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;
    
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

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

    /**
     * @brief It calculates the constitutive matrix C
     * @param C The constitutive matrix
     * @param rValues Parameters of the constitutive law
     */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& C,
        ConstitutiveLaw::Parameters& rValues
        ) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(
        const ConstitutiveLaw::StrainVectorType& rStrainVector,
        ConstitutiveLaw::StressVectorType& rStressVector,
        ConstitutiveLaw::Parameters& rValues
        ) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallStrainJ2Elasticity3D)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallStrainJ2Elasticity3D)
    }


}; // Class SmallStrainJ2ElasticityPlaneStrain
}  // namespace Kratos.
