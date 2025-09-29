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
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Alfredo Huespe
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "small_strain_phase_trans_thermal_plasticity_3d.h"
#include "constitutive_laws_small_strain_application_variables.h"
#include "structural_analysis_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainPhaseTransThermalPlasticity3D::SmallStrainPhaseTransThermalPlasticity3D()
    : SmallStrainJ2ThermalPlasticity3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainPhaseTransThermalPlasticity3D::SmallStrainPhaseTransThermalPlasticity3D(const SmallStrainPhaseTransThermalPlasticity3D &rOther)
    : SmallStrainJ2ThermalPlasticity3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainPhaseTransThermalPlasticity3D::Clone() const
{
    return Kratos::make_shared<SmallStrainPhaseTransThermalPlasticity3D>(SmallStrainPhaseTransThermalPlasticity3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainPhaseTransThermalPlasticity3D::~SmallStrainPhaseTransThermalPlasticity3D()
{
}

//************************************************************************************
//************************************************************************************




//************************************************************************************
//************************************************************************************

/*++
subroutine to return the phase transformation strain component increment, 
mfinc[0] = austenite increment,
mfinc[1,2,3] = eutectoid phase incrememnt,
mfinc[4,5] = bainite increment,
mfinc[6] = martensite increment.
--*/
double SmallStrainPhaseTransThermalPlasticity3D::GetPhaseTransStrainInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc)
{
	const auto& rMaterialProperties = rValues.GetMaterialProperties();
    const auto& rElementGeometry = rValues.GetElementGeometry();
    const auto& rCurrentProcessInfo = rValues.GetProcessInfo();
    
    double dret = 0;

    if (mfinc[0] > 0 && this->mpPT[0] != nullptr) {
        dret = this->mpPT[0]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[0];
    }
    else {
        if (this->mpPT[1] != nullptr) dret = this->mpPT[1]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[1] + mfinc[2] + mfinc[3]);
        if (this->mpPT[2] != nullptr) dret += this->mpPT[2]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[4] + mfinc[5]);
        if (this->mpPT[3] != nullptr) dret += this->mpPT[3]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[6];
    }

	/* return the data */
	return dret;
}

/*++
subroutine to return the creep strain increment
--*/
int SmallStrainPhaseTransThermalPlasticity3D::GetCreepStrainInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc)
{
	const auto& rMaterialProperties = rValues.GetMaterialProperties();
    const auto& rElementGeometry = rValues.GetElementGeometry();
    const auto& rCurrentProcessInfo = rValues.GetProcessInfo();
    
    double dret = 0;

    if (mfinc[0] > 0 && this->mpPT[0] != nullptr) {
        dret = this->mpPT[0]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[0];
    }
    else {
        if (this->mpPT[1] != nullptr) dret = this->mpPT[1]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[1] + mfinc[2] + mfinc[3]);
        if (this->mpPT[2] != nullptr) dret += this->mpPT[2]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[4] + mfinc[5]);
        if (this->mpPT[3] != nullptr) dret += this->mpPT[3]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[6];
    }

	/* return the data */
	return dret;
}

/*++
subroutine to return the integral mean value of creep strain increment
--*/
int SmallStrainPhaseTransThermalPlasticity3D::GetCreepStrainIncMean(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc)
{
	const auto& rMaterialProperties = rValues.GetMaterialProperties();
    const auto& rElementGeometry = rValues.GetElementGeometry();
    const auto& rCurrentProcessInfo = rValues.GetProcessInfo();
    
    double dret = 0;

    if (mfinc[0] > 0 && this->mpPT[0] != nullptr) {
        dret = this->mpPT[0]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[0];
    }
    else {
        if (this->mpPT[1] != nullptr) dret = this->mpPT[1]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[1] + mfinc[2] + mfinc[3]);
        if (this->mpPT[2] != nullptr) dret += this->mpPT[2]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[4] + mfinc[5]);
        if (this->mpPT[3] != nullptr) dret += this->mpPT[3]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[6];
    }

	/* return the data */
	return dret;
}

/*++
subroutine to return the trip increment,
mf[0], mfinc[0] = austenite and its increment,
mf[1], mfinc[1] = eutectoid phases and their incrememnt,
mf[2], mfinc[2] = bainite and their increment,
mf[3], mfinc[3] = martensite and its increment.
--*/
int SmallStrainPhaseTransThermalPlasticity3D::GetTripInc(Kratos::BasePhaseTransitionModel::Parameters &rValues, double *mfinc)
{
	const auto& rMaterialProperties = rValues.GetMaterialProperties();
    const auto& rElementGeometry = rValues.GetElementGeometry();
    const auto& rCurrentProcessInfo = rValues.GetProcessInfo();
    
    double dret = 0;

    if (mfinc[0] > 0 && this->mpPT[0] != nullptr) {
        dret = this->mpPT[0]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[0];
    }
    else {
        if (this->mpPT[1] != nullptr) dret = this->mpPT[1]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[1] + mfinc[2] + mfinc[3]);
        if (this->mpPT[2] != nullptr) dret += this->mpPT[2]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * (mfinc[4] + mfinc[5]);
        if (this->mpPT[3] != nullptr) dret += this->mpPT[3]->GetPhaseTransStrainMean(rMaterialProperties, rElementGeometry, rCurrentProcessInfo) * mfinc[6];
    }

	/* return the data */
	return dret;
}

void SmallStrainPhaseTransThermalPlasticity3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->GetStrainSize();
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************

int SmallStrainPhaseTransThermalPlasticity3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT));

    return 0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainPhaseTransThermalPlasticity3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.save("mThermalStrainNormal", mThermalStrainNormal);
    rSerializer.save("mThermalStrain", mThermalStrain);
}

//************************************************************************************
//************************************************************************************

void SmallStrainPhaseTransThermalPlasticity3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.load("mThermalStrainNormal", mThermalStrainNormal);
    rSerializer.load("mThermalStrain", mThermalStrain);
}

} /* namespace Kratos.*/
