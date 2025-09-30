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
#include "../../elasticity/small_strain_j2_elasticity_3d.h"
#include "phase_transformation/phase_transformation_system.h"

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
 * @class SmallStrainJ2ThermalElasticity3D
 * @ingroup constitutive_laws_small_strain_application
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 */
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) SmallStrainJ2ThermalElasticity3D
    : public SmallStrainJ2Elasticity3D
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = SmallStrainJ2Elasticity3D;

    /// Counted pointer of LinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainJ2ThermalElasticity3D);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainJ2ThermalElasticity3D() 
    {}

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SmallStrainJ2ThermalElasticity3D>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~SmallStrainJ2ThermalElasticity3D() override {}

    /**
    * Copy constructor.
    */
    SmallStrainJ2ThermalElasticity3D(const SmallStrainJ2ThermalElasticity3D &rOther)
        : BaseType(rOther),
        mThermalStrain(rOther.mThermalStrain)
        mThermalStrainNormal(rOther.mThermalStrainNormal)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function set phase transformation system pointer
     */
    virtual void SetPhaseTransSystem(PhaseTransformationSystem* pSystem) { mpPhaseTransSystem = pSystem }

    /**
     * @brief This function get phase transformation system pointer
     */
    const PhaseTransformationSystem* GetPhaseTransSystem() const { return mpPhaseTransSystem; }

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties &rMaterialProperties,
        const GeometryType &rElementGeometry,
        const Vector &rShapeFunctionsValues) override;

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;


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
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;


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
    virtual void CalculateThermalStrainIncrement(Vector& rThermalStrainIncrementVector, ConstitutiveLaw::Parameters &rValues);


    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
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


    virtual double& GetThermalStrainNormal()    { return mThermalStrainNormal; }
    virtual Vector& GetThermalStrain()          { return mThermalStrain; }
    virtual void SetThermalStrainNormal(const double ToThermalStrainNormal) {mThermalStrainNormal = ToThermalStrainNormal;}
    virtual void SetThermalStrain(const array_1d<double, 6>& rThermalStrain) {mThermalStrain = rThermalStrain;}

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

    Vector mThermalStrain;
    double mThermalStrainNormal;
    PhaseTransformationSystem* mpPhaseTransSystem

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

    /**
     * @brief Retrieve the reference temperature
     * @return The reference temperature
     */

    /**
     * @brief Sets the reference temperature
     * @param ToRefTemperature The reference temperature
     */

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainJ2Elasticity3D)
        rSerializer.save("ThermalStrain", mThermalStrain);
        rSerializer.save("ThermalStrainNormal", mThermalStrainNormal);
        rSerializer.save("PhaseTransSystem", PhaseTransSystem);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainJ2Elasticity3D)
        rSerializer.load("ThermalStrain", mThermalStrain);
        rSerializer.load("ThermalStrainNormal", mThermalStrainNormal);
        rSerializer.load("PhaseTransSystem", PhaseTransSystem);
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
