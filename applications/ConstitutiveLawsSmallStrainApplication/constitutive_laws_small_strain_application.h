// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//

#pragma once


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Plasticity Constitutive laws
#include "custom_constitutive/small_strains/plasticity/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_plane_strain_2d.h"

// Integrators
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_kinematic_plasticity.h"


// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"

// Thermal yield surfaces
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"

// Rules of mixtures

// Thermal CL
#include "custom_constitutive/thermal/small_strains/elastic/thermal_elastic_isotropic_3d.h"
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_strain.h"
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_stress.h"
#include "custom_constitutive/thermal/small_strains/plastic/thermal_plastic_isotropic_3d.h"
#include "custom_constitutive/thermal/small_strains/plastic/thermal_plastic_linear_plane_strain.h"
#include "custom_constitutive/thermal/small_strains/plastic/thermal_plastic_linear_plane_stress.h"
#include "custom_constitutive/thermal/small_strains/plastic/generic_small_strain_thermal_isotropic_plasticity.h"
#include "custom_constitutive/thermal/small_strains/plastic/generic_small_strain_thermal_isotropic_plasticity_plane_stress.h"
#include "custom_constitutive/thermal/small_strains/plastic/small_strain_j2_thermal_plasticity_3d.h"
#include "custom_constitutive/thermal/small_strains/plastic/small_strain_j2_thermal_plasticity_plane_strain_2d.h"

namespace Kratos {

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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CONSTITUTIVE_LAWS_SMALL_STRAIN_APPLICATION) KratosConstitutiveLawsSmallStrainApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosConstitutiveLawsSmallStrainApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosConstitutiveLawsSmallStrainApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConstitutiveLawsSmallStrainApplication();

    /// Destructor.
    ~KratosConstitutiveLawsSmallStrainApplication() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosConstitutiveLawsSmallStrainApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
          KRATOS_WATCH("in my application");
          KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

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

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{


    // Plasticity laws
    const SmallStrainIsotropicPlasticityFactory mSmallStrainIsotropicPlasticityFactory;

    //Damage laws

    /// Plasticity
    /* Small strain */
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesVonMises;
    const SmallStrainJ2Plasticity3D mSmallStrainJ2Plasticity3D;
    const SmallStrainJ2PlasticityPlaneStrain2D mSmallStrainJ2PlasticityPlaneStrain2D;

    // Plastic Damage Model


    /* Finite strain Isotropic plasticity*/


    /* Finite strain Kinematic plasticity*/


    /// Damage
    /* Small strain 3D */

    /* Small strain 2D */


    /* Small strain  plane stress 2D */


    // HCF (High Cycle Fatigue)



    // d+d- laws (3D)


    // d+d- laws (2D)


    // Orthotropic Damage


    // Rules of mixtures


    // Anisotropic law


    // Thermal CL
    const ThermalElasticIsotropic3D mThermalElasticIsotropic3D;
    const ThermalLinearPlaneStrain mThermalLinearPlaneStrain;
    const ThermalLinearPlaneStress mThermalLinearPlaneStress;
    const ThermalPlasticIsotropic3D mThermalPlasticIsotropic3D;
    const ThermalPlasticLinearPlaneStrain mThermalPlasticLinearPlaneStrain;
    const ThermalPlasticLinearPlaneStress mThermalPlasticLinearPlaneStress;
    const GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicPlasticity3DVonMisesVonMises;
    const GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicPlasticityPlaneStressVonMisesVonMises;

    const SmallStrainJ2ThermalPlasticity3D mSmallStrainJ2ThermalPlasticity3D;
    const SmallStrainJ2ThermalPlasticityPlaneStrain2D mSmallStrainJ2ThermalPlasticityPlaneStrain2D;

    /// Thermal Damage
    /* Small strain 3D */


    /* Small strain 2D */

    // Thermal plane stress damage


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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosConstitutiveLawsSmallStrainApplication& operator=(KratosConstitutiveLawsSmallStrainApplication const& rOther);

    /// Copy constructor.
    KratosConstitutiveLawsSmallStrainApplication(KratosConstitutiveLawsSmallStrainApplication const& rOther);


    ///@}

}; // Class KratosConstitutiveLawsSmallStrainApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.
