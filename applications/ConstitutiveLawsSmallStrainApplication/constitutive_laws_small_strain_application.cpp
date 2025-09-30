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


// System includes


// External includes


// Project includes
#include "constitutive_laws_small_strain_application.h"
#include "constitutive_laws_small_strain_application_variables.h"


namespace Kratos {

KratosConstitutiveLawsSmallStrainApplication::KratosConstitutiveLawsSmallStrainApplication():
    KratosApplication("ConstitutiveLawsSmallStrainApplication")
    {}

void KratosConstitutiveLawsSmallStrainApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosConstitutiveLawsSmallStrainApplication..." << std::endl;

    /// Elasticity
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2Elasticity3D", mSmallStrainJ2Elasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ElasticityPlaneStrain", mSmallStrainJ2ElasticityPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ElasticityPlaneStress", mSmallStrainJ2ElasticityPlaneStress);

    /// Plasticity

    /* Small strain */
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2Plasticity3D", mSmallStrainJ2Plasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2PlasticityPlaneStrain", mSmallStrainJ2PlasticityPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2PlasticityPlaneStress", mSmallStrainJ2PlasticityPlaneStress);

    // Thermal CL
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ThermalElasticity3D", mSmallStrainJ2ThermalElasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ThermalElasticityPlaneStrain", mSmallStrainJ2ThermalElasticityPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ThermalElasticityPlaneStress", mSmallStrainJ2ThermalElasticityPlaneStress);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ThermalPlasticity3D", mSmallStrainJ2ThermalPlasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ThermalPlasticityPlaneStrain", mSmallStrainJ2ThermalPlasticityPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2ThermalPlasticityPlaneStress", mSmallStrainJ2ThermalPlasticityPlaneStress);


    // Constitutive laws variables

    KRATOS_REGISTER_VARIABLE(MAX_NUMBER_NL_CL_ITERATIONS)

    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_TENSOR)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR)

    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS)
    
    KRATOS_REGISTER_VARIABLE(YIELD_STRENGTH)

    KRATOS_REGISTER_VARIABLE(THERMAL_STRAIN_VECTOR)
    KRATOS_REGISTER_VARIABLE(THERMAL_STRAIN_NORMAL)

}
}  // namespace Kratos.
