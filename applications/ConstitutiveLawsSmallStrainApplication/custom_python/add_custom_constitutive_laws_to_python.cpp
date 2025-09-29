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
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Elastic laws


// Plastic laws
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_plane_strain_2d.h"

//Damage laws

//Viscosities

// Integrators
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
// Thermal yield surfaces
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"


// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"

// Rules of mixtures


// Thermal CL's
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

#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"

namespace Kratos::Python {

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;



    py::class_< SmallStrainIsotropicPlasticityFactory, typename SmallStrainIsotropicPlasticityFactory::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicPlasticityFactory").def(py::init<>())
    ;


    // Custom Constitutive Laws Registration
    // Isotropic Plasticity
    /* Small strain */

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises").def(py::init<>());


    // Isotropic Plasticity 3D
    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< SmallStrainJ2Plasticity3D, typename SmallStrainJ2Plasticity3D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2Plasticity3DLaw").def(py::init<>());

    py::class_< SmallStrainJ2PlasticityPlaneStrain2D, typename SmallStrainJ2PlasticityPlaneStrain2D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2PlasticityPlaneStrain2D").def(py::init<>());


    // Plastic Damage Model


    // Kinematic Plasticity
    /* Small strain */


    // Kinematic isotropic 3D plasticity

    // HCF

    // kinematic plasticity

    /* Finite strain */
    // Isotropic


    /* Finite strain */
    // Kinematic


    // Damage
    /* Small strain */

    // damage 3D

    // damage plane strain

    // damage plane stress



    // Thermal CL's
    py::class_< ThermalElasticIsotropic3D, typename ThermalElasticIsotropic3D::Pointer, ConstitutiveLaw >
    (m,"ThermalElasticIsotropic3D").def(py::init<>());

    py::class_< ThermalLinearPlaneStrain, typename ThermalLinearPlaneStrain::Pointer, ConstitutiveLaw >
    (m,"ThermalLinearPlaneStrain").def(py::init<>());

    py::class_< ThermalLinearPlaneStress, typename ThermalLinearPlaneStress::Pointer, ConstitutiveLaw >
    (m,"ThermalLinearPlaneStress").def(py::init<>());
    py::class_< ThermalPlasticIsotropic3D, typename ThermalPlasticIsotropic3D::Pointer, ConstitutiveLaw >
    (m,"ThermalPlasticIsotropic3D").def(py::init<>());

    py::class_< ThermalPlasticLinearPlaneStrain, typename ThermalPlasticLinearPlaneStrain::Pointer, ConstitutiveLaw >
    (m,"ThermalPlasticLinearPlaneStrain").def(py::init<>());

    py::class_< ThermalPlasticLinearPlaneStress, typename ThermalPlasticLinearPlaneStress::Pointer, ConstitutiveLaw >
    (m,"ThermalPlasticLinearPlaneStress").def(py::init<>());
    
    py::class_<  GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicPlasticity3DVonMises").def(py::init<>());
    
    py::class_<  GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicPlasticityPlaneStrainVonMises").def(py::init<>());

    py::class_< SmallStrainJ2ThermalPlasticity3D, typename SmallStrainJ2ThermalPlasticity3D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalPlasticity3DLaw").def(py::init<>());

    py::class_< SmallStrainJ2ThermalPlasticityPlaneStrain2D, typename SmallStrainJ2ThermalPlasticityPlaneStrain2D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalPlasticityPlaneStrain2DLaw").def(py::init<>());


    // damage 3D

    // damage plane strain

    // Plane stress thermal damage


}

}  // namespace Kratos::Python.
