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

namespace Kratos::Python {

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;



    // Custom Constitutive Laws Small strain Registration

    // Isotropic Elasticity
    py::class_< SmallStrainJ2Elasticity3D, typename SmallStrainJ2Elasticity3D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2Elasticity3D").def(py::init<>());

    py::class_< SmallStrainJ2ElasticityPlaneStrain, typename SmallStrainJ2ElasticityPlaneStrain::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2ElasticityPlaneStrain").def(py::init<>());

    py::class_< SmallStrainJ2ElasticityPlaneStress, typename SmallStrainJ2ElasticityPlaneStress::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2ElasticityPlaneStress").def(py::init<>());


    // Isotropic Plasticity 3D
    py::class_< SmallStrainJ2Plasticity3D, typename SmallStrainJ2Plasticity3D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2Plasticity3D").def(py::init<>());

    py::class_< SmallStrainJ2PlasticityPlaneStrain, typename SmallStrainJ2PlasticityPlaneStrain::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2PlasticityPlaneStrain").def(py::init<>());

    py::class_< SmallStrainJ2PlasticityPlaneStress, typename SmallStrainJ2PlasticityPlaneStress::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2PlasticityPlaneStress").def(py::init<>());


    // Thermal CL's
    py::class_< SmallStrainJ2ThermalElasticity3D, typename SmallStrainJ2ThermalElasticity3D::Pointer, ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalElasticity3D").def(py::init<>());

    py::class_< SmallStrainJ2ThermalElasticityPlaneStrain, typename SmallStrainJ2ThermalElasticityPlaneStrain::Pointer, ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalElasticityPlaneStrain").def(py::init<>());

    py::class_< SmallStrainJ2ThermalElasticityPlaneStress, typename SmallStrainJ2ThermalElasticityPlaneStress::Pointer, ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalElasticityPlaneStress").def(py::init<>());

    py::class_< SmallStrainJ2ThermalPlasticity3D, typename SmallStrainJ2ThermalPlasticity3D::Pointer, ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalPlasticity3D").def(py::init<>());

    py::class_< SmallStrainJ2ThermalPlasticityPlaneStrain, typename SmallStrainJ2ThermalPlasticityPlaneStrain::Pointer, ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalPlasticityPlaneStrain").def(py::init<>());

    py::class_< SmallStrainJ2ThermalPlasticityPlaneStress, typename SmallStrainJ2ThermalPlasticityPlaneStress::Pointer, ConstitutiveLaw >
    (m,"SmallStrainJ2ThermalPlasticityPlaneStress").def(py::init<>());
    


}

}  // namespace Kratos::Python.
