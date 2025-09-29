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
//                  whf

#pragma once

// System includes

// External includes
// Project includes
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"
#include "custom_utilities/advanced_constitutive_law_small_strain_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    typedef std::size_t SizeType;

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
 * @class GenericSmallStrainThermalIsotropicPlasticity
 * @ingroup ConstitutiveLawsApp
 * @brief This class derives from the Isotropic plasticity CL and adds thermal effects (material properties affectation and internal variables)
 * @details This class considers a constitutive law integrator as an intermediate utility to compute the plasticity. 3D and plane strain
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author whf
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) GenericSmallStrainThermalIsotropicPlasticity
    : public GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TConstLawIntegratorType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TConstLawIntegratorType::VoigtSize;

    /// Definition of the base class
    using BaseType = GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>;

    /// Bounded vector for stresses/strains
    using BoundedArrayType = array_1d<double, VoigtSize>;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainThermalIsotropicPlasticity);

    /// The geometry definition
    using GeometryType = Geometry<Node>;

    /// Advanced and basic contitutive laws utilities for the corresponding Voigt size
    using CLutils    = ConstitutiveLawUtilities<VoigtSize>;
    using AdvCLutils = AdvancedConstitutiveLawSmallStrainUtilities<VoigtSize>;

    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();
    static constexpr double threshold_tolerance = 1.0e-5;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainThermalIsotropicPlasticity()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainThermalIsotropicPlasticity<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainThermalIsotropicPlasticity(const GenericSmallStrainThermalIsotropicPlasticity &rOther)
        : BaseType(rOther),
          mReferenceTemperature(rOther.mReferenceTemperature),
          mThermalStrain(rOther.mThermalStrain),
          mThermalStrainNormal(rOther.mThermalStrainNormal)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainThermalIsotropicPlasticity() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Calculate the thermal strain vector increment
     * @param rThermalStrainVectorIncrement The thermal strain vector increment
     * @param rValues The constitutive law parameters
     * @param rReferenceTemperature The reference temperature
     */
    virtual void CalculateThermalStrainIncrement(Vector& rThermalStrainIncrementVector, ConstitutiveLaw::Parameters &rValues);

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double> &rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector> &rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double> &rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector> &rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double> &rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector> &rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;


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
    ///@name Protected  Access
    ///@{

    /**
     * @brief Retrieve the reference temperature
     * @return The reference temperature
     */
    double& GetReferenceTemperature()
    {
        return mReferenceTemperature;
    }

    /**
     * @brief Sets the reference temperature
     * @param ToRefTemperature The reference temperature
     */
    void SetReferenceTemperature(const double ToRefTemperature)
    {
        mReferenceTemperature = ToRefTemperature;
    }

    double& GetThermalStrainNormal()
    {
        return mThermalStrainNormal;
    }

    /**
     * @brief Sets the thermal strain normal
     * @param ToThermalStrainNormal The thermal strain normal
     */
    void SetThermalStrainNormal(const double ToThermalStrainNormal)
    {
        mThermalStrainNormal = ToThermalStrainNormal;
    }

    Vector& GetThermalStrain()
    {
        return mThermalStrain;
    }

    void SetThermalStrain(const array_1d<double, VoigtSize>& rThermalStrain)
    {
        mThermalStrain = rThermalStrain;
    }
    
    void CalculateTangentTensor(
        ConstitutiveLaw::Parameters &rValues,
        const Vector& rPlasticStrain);

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

    double mReferenceTemperature = 0.0;
    
    double mPlasticDissipation = 0.0;
    double mThreshold = 0.0;
    Vector mPlasticStrain = ZeroVector(VoigtSize);
    Vector mThermalStrain = ZeroVector(VoigtSize);
    double mThermalStrainNormal = 0.0;

    ///@}
    ///@name Private Operators
    ///@{
    
   /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>)
        rSerializer.save("ReferenceTemperature", mReferenceTemperature);
        rSerializer.save("ThermalStrain", mThermalStrain);
        rSerializer.save("ThermalStrainNormal", mThermalStrainNormal);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>)
        rSerializer.load("ReferenceTemperature", mReferenceTemperature);
        rSerializer.load("ThermalStrain", mThermalStrain);
        rSerializer.load("ThermalStrainNormal", mThermalStrainNormal);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
