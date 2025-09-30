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
//  Main authors:    whf
//

#pragma once


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Elasticity Constitutive laws
#include "custom_constitutive/elasticity/small_strain_j2_elasticity_3d"
#include "custom_constitutive/elasticity/small_strain_j2_elasticity_plane_strain"
#include "custom_constitutive/elasticity/small_strain_j2_elasticity_plane_stress"

// Plasticity Constitutive laws
#include "custom_constitutive/plasticity/small_strain_j2_plasticity_factory.h"
#include "custom_constitutive/plasticity/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/plasticity/small_strain_j2_plasticity_plane_strain.h"
#include "custom_constitutive/plasticity/small_strain_j2_plasticity_plane_stress.h"

// Rules of mixtures

// Thermal CL
#include "custom_constitutive/thermal/elastic/small_strain_j2_thermal_elasticity_3d.h"
#include "custom_constitutive/thermal/elastic/small_strain_j2_thermal_elasticity_plane_strain.h"
#include "custom_constitutive/thermal/elastic/small_strain_j2_thermal_elasticity_plane_stress.h"
#include "custom_constitutive/thermal/plastic/small_strain_j2_thermal_plasticity_3d.h"
#include "custom_constitutive/thermal/plastic/small_strain_j2_thermal_plasticity_plane_strain.h"
#include "custom_constitutive/thermal/plastic/small_strain_j2_thermal_plasticity_plane_stress.h"


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


    // Elasticity laws
    const SmallStrainJ2Elasticity3D mSmallStrainJ2Elasticity3D;
    const SmallStrainJ2ElasticityPlaneStrain mSmallStrainJ2ElasticityPlaneStrain;
    const SmallStrainJ2ElasticityPlaneStress mSmallStrainJ2ElasticityPlaneStress;
    
    // Plasticity laws
    const SmallStrainJ2Plasticity3D mSmallStrainJ2Plasticity3D;
    const SmallStrainJ2PlasticityPlaneStrain mSmallStrainJ2PlasticityPlaneStrain;
    const SmallStrainJ2PlasticityPlaneStress mSmallStrainJ2PlasticityPlaneStress;

    // Thermal CL
    const SmallStrainJ2ThermalElasticity3D mSmallStrainJ2ThermalElasticity3D;
    const SmallStrainJ2ThermalElasticityPlaneStrain mSmallStrainJ2ThermalElasticityPlaneStrain;
    const SmallStrainJ2ThermalElasticityPlaneStress mSmallStrainJ2ThermalElasticityPlaneStress;
    const SmallStrainJ2ThermalPlasticity3D mSmallStrainJ2ThermalPlasticity3D;
    const SmallStrainJ2ThermalPlasticityPlaneStrain mSmallStrainJ2ThermalPlasticityPlaneStrain;
    const SmallStrainJ2ThermalPlasticityPlaneStress mSmallStrainJ2ThermalPlasticityPlaneStress;


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
