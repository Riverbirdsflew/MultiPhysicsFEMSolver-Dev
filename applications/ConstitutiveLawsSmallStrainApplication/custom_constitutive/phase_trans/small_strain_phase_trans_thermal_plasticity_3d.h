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
//  Main authors:    whf
//                   

#pragma once

// System includes

// External includes

// Project includes
#include "../thermal/small_strains/plastic/small_strain_j2_thermal_plasticity_3d.h"
#include "custom_utilities/advanced_constitutive_law_small_strain_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "PhaseTransitionsApplication/phase_transition/phase_transition_model/base_phase_transition_model.h"

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
 * @class J2ThermalPlasticity3D
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a Simo J2 thermal plasticity constitutive law in 3D
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - YIELD_STRESS
 * @warning Valid for small strains, linear hexahedra
 * @author whf
 */
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) SmallStrainPhaseTransThermalPlasticity3D
    : public SmallStrainJ2ThermalPlasticity3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    using CLutils    = ConstitutiveLawUtilities<6>;
    using AdvCLutils = AdvancedConstitutiveLawSmallStrainUtilities<6>;

    // Counted pointer of SmallStrainPhaseTransThermalPlasticity3D
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainPhaseTransThermalPlasticity3D);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SmallStrainPhaseTransThermalPlasticity3D();

    /**
     * @brief Copy constructor.
     */
    SmallStrainPhaseTransThermalPlasticity3D(const SmallStrainPhaseTransThermalPlasticity3D& rOther);

    /**
     * @brief Destructor.
     */
    ~SmallStrainPhaseTransThermalPlasticity3D() override;

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
     * @brief dimension of the constitutive law
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return 6;
    };

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rThisVariable The variable to be returned
     * @param rValue New value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable The variable to be returned
     * @param rValue New value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector>& rThisVariable,
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
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

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
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

	/* return the phase transformation strain component increment */
	virtual double GetPhaseTransStrainInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the creep cofficient data */
	virtual int GetCreepParameters(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the creep strain increment */
	virtual int GetCreepStrainInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the integral mean value of creep strain increment */
	virtual int GetCreepStrainIncMean(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the trip increment */
	virtual int GetTripInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);

	/* return the phase transformation latent */
	//virtual double GetLatent(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);

	/* phase transformation part */
	/* return the austenite mass fraction increment */
	virtual int GetAusteniteMassFrcInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the eutectoid phase mass fraction increment */
	virtual int GetEutectoidPhaseMassFrcInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the bainite mass fraction increment */
	virtual int GetBainiteMassFrcInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* return the martensite mass fraction increment */
	virtual int GetMartensiteMassFrcInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
	/* update the incubation, carbon content and mass fraction of all microstructures */
	virtual int UpdateMassFraction(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc);
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
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SmallStrainPhaseTransThermalPlasticity3D";
        return buffer.str();
    }
    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "J2 Thermal Plasticity 3D constitutive law\n";
    };

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    BasePhaseTransitionModel::Pointer mpPT[4]; // Pointer to the phase transition model
    double mTransitionStrain;
    Vector mTripStrain;

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
    using ConstitutiveLaw::CalculateStressResponse;
    virtual void CalculateStressResponse(ConstitutiveLaw::Parameters& rValues,
                                 Vector& rPlasticStrain,
                                 double& rAccumulatedPlasticStrain,
                                 Vector& rThermalStrain);

    /**
     * @brief This method computes the yield function
     * @param NormDeviationStress The norm of the deviation stress
     * @param rMaterialProperties The properties of the current material considered
     * @return The trial yield function (after update)
     */
    double YieldFunction(
        const double NormDeviationStress,
        ConstitutiveLaw::Parameters& rValues,
        const double AccumulatedPlasticStrain,
        const double Dgamma = 0.0
        );

    /**
     * @brief This method computes the hardening modulus
     * @param NormDeviationStress The norm of the deviation stress
     * @param rMaterialProperties The properties of the current material considered
     * @return The trial yield function (after update)
     */
    double GetHardeningModulus(
        ConstitutiveLaw::Parameters& rValues,
        const double AccumulatedPlasticStrain
        );

    /**
     * @brief This method computes the increment of Gamma
     * @param NormStressTrial The norm of the stress trial
     * @param rMaterialProperties The properties of the material
     * @return The increment of Gamma computed
     */
    double GetDgamma(const double NormStressTrial, ConstitutiveLaw::Parameters& rValues,
                                     const double AccumulatedPlasticStrainOld);


 /**
     * @brief This method gets the saturation hardening parameter
     * @param rMaterialProperties The properties of the material
     * @return The saturation hardening parameter
     */
    double GetYieldStress(ConstitutiveLaw::Parameters& rValues, const double AccumulatedPlasticStrain);


    /**
     * @brief This method computes the constitutive tensor
     * @param DeltaGamma The increment on the Gamma parameter
     * @param NormStressTrial The norm of the stress trial
     * @param rYFNormalVector The yield function normal vector
     * @param rMaterialProperties The properties of the material
     * @param rTMatrix The elastic tensor/matrix to be computed
     */
    virtual void CalculateTangentMatrix(const double Dgamma, const double NormTrialStressDev,
                                        const Vector &rYFNormalVector,
                                        ConstitutiveLaw::Parameters &rValues,
                                        const double AccumulatedPlasticStrain, Matrix &rTMatrix);

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticMatrix The elastic tensor/matrix to be computed
     * @param rMaterialProperties The properties of the material
     */
    virtual void CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues, Matrix &rElasticMatrix);

    ///@}
    ///@name Protected  Access
    ///@{

    virtual double& GetThermalStrainNormal()    { return mThermalStrainNormal; }
    virtual Vector& GetThermalStrain()          { return mThermalStrain; }
    virtual void SetThermalStrainNormal(const double ToThermalStrainNormal) {mThermalStrainNormal = ToThermalStrainNormal;}
    virtual void SetThermalStrain(const array_1d<double, 6>& rThermalStrain) {mThermalStrain = rThermalStrain;}


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

    ///@}
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

}; // Class SmallStrainPhaseTransThermalPlasticity3D
} // namespace Kratos.
