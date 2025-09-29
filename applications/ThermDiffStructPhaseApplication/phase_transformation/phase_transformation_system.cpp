//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    WHF
//

// System includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <chrono>

// External includes

// Project includes
#include "phase_transformation_system.h"
#include "../therm_diff_struct_phase_application_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

PhaseTransformationSystem::PhaseTransformationSystem()
{
    SizeType i;

    this->mstat = Status::UnSet;

    for(i=0; i<NumPhases; ++i) {
        this->mpPhases[i] = nullptr;
    }
    for(i=0; i<NumTransformations; ++i) {
        this->mpPhaseTrans[i] = nullptr;
    } 
}

PhaseTransformationSystem::~PhaseTransformationSystem()
{
    SizeType i;
    for(i=0; i<NumPhases; ++i) 
    {
        if(this->mpPhases[i] != nullptr) 
        {
            this->mpPhases[i].reset();
            this->mpPhases[i] = nullptr;
        }
    }
    for(i=0; i<NumTransformations; ++i) {
        if(this->mpPhases[i] != nullptr) 
        {
            this->mpPhases[i].reset();
            this->mpPhases[i] = nullptr;
        }
    } 
    this->mstat = Status::UnSet;
}


/**
 * 返回内部对象
 */
const PhaseTransformationSystem::Status PhaseTransformationSystem::Stat()
{
	/* return the flag of data status */
	return this->mstat;
}

const Phase::Pointer PhaseTransformationSystem::GetPhaseObjectPtr(unsigned int id)
{
    if (id < 0 || id >= NumPhases) {
        KRATOS_ERROR << "Invalid PhaseIndex: " << id << std::endl;
    }
    return this->mpPhases[id];
}

const PhaseTransformationLaw::Pointer PhaseTransformationSystem::GetPhaseTransObjectPtr(unsigned int id)
{
    if (id < 0 || id >= NumTransformations) {
        KRATOS_ERROR << "Invalid PhaseIndex: " << id << std::endl;
    }
    return this->mpPhaseTrans[id];
}

//************************************************************************************
//************************************************************************************

int PhaseTransformationSystem::Initialize(const Properties& rMaterialProperties,
                                            PhaseTransformationType rTransType,
                                            PhaseType rPhaseType)
{
    int iret = 0x0;

    KRATOS_TRY;

    if(rMaterialProperties.Has(PHASE_TRANSFORMATION_LAW_AUSTENIZATION))
    {
        auto& r_var = KratosComponents<Variable<PhaseTransformationLaw::Pointer>>::Get("PHASE_TRANSFORMATION_LAW_AUSTENIZATION");
        this->mpPhaseTrans[0] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhaseTrans[0] == nullptr) << "PhaseTransformationSystem:: PHASE_TRANSFORMATION_LAW_AUSTENIZATION is not set!" << std::endl;
        this->mpPhaseTrans[0]->Initialize(rMaterialProperties, rTransType);
    }
    if(rMaterialProperties.Has(PHASE_TRANSFORMATION_LAW_EUTECTOID_DECOMPOSITION))
    {
        auto& r_var = KratosComponents<Variable<PhaseTransformationLaw::Pointer>>::Get("PHASE_TRANSFORMATION_LAW_EUTECTOID_DECOMPOSITION");
        this->mpPhaseTrans[1] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhaseTrans[1] == nullptr) << "PhaseTransformationSystem:: PHASE_TRANSFORMATION_LAW_EUTECTOID_DECOMPOSITION is not set!" << std::endl;
        this->mpPhaseTrans[1]->Initialize(rMaterialProperties, rTransType);
    }
    if(rMaterialProperties.Has(PHASE_TRANSFORMATION_LAW_BAINITIC_TRANSFORMATION))
    {
        auto& r_var = KratosComponents<Variable<PhaseTransformationLaw::Pointer>>::Get("PHASE_TRANSFORMATION_LAW_BAINITIC_TRANSFORMATION");
        this->mpPhaseTrans[2] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhaseTrans[2] == nullptr) << "PhaseTransformationSystem:: PHASE_TRANSFORMATION_LAW_BAINITIC_TRANSFORMATION is not set!" << std::endl;
        this->mpPhaseTrans[2]->Initialize(rMaterialProperties, rTransType);
    }
    if(rMaterialProperties.Has(PHASE_TRANSFORMATION_LAW_MARTENSITIC_TRANSFORMATION))
    {
        auto& r_var = KratosComponents<Variable<PhaseTransformationLaw::Pointer>>::Get("PHASE_TRANSFORMATION_LAW_MARTENSITIC_TRANSFORMATION");
        this->mpPhaseTrans[3] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhaseTrans[3] == nullptr) << "PhaseTransformationSystem:: PHASE_TRANSFORMATION_LAW_MARTENSITIC_TRANSFORMATION is not set!" << std::endl;
        this->mpPhaseTrans[3]->Initialize(rMaterialProperties, rTransType);
    }
    if(rMaterialProperties.Has(PHASE_AUSTENITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_AUSTENITE");
        this->mpPhases[0] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[0] == nullptr) << "PhaseTransformationSystem:: AUSTENITE is not set!" << std::endl;
        this->mpPhases[0]->Initialize(rMaterialProperties, rPhaseType);
    }
    if(rMaterialProperties.Has(PHASE_FERRITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_FERRITE");
        this->mpPhases[1] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[1] == nullptr) << "PhaseTransformationSystem:: FERRITE is not set!" << std::endl;
        this->mpPhases[1]->Initialize(rMaterialProperties, rPhaseType);
    }
    if(rMaterialProperties.Has(PHASE_CEMENTITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_CEMENTITE");
        this->mpPhases[2] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[2] == nullptr) << "PhaseTransformationSystem:: CEMENTITE is not set!" << std::endl;
        this->mpPhases[2]->Initialize(rMaterialProperties, rPhaseType);
    }
    if(rMaterialProperties.Has(PHASE_PEARLITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_PEARLITE");
        this->mpPhases[3] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[3] == nullptr) << "PhaseTransformationSystem:: PEARLITE is not set!" << std::endl;
        this->mpPhases[3]->Initialize(rMaterialProperties, rPhaseType);
    }
    if(rMaterialProperties.Has(PHASE_UPPER_BAINITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_UPPER_BAINITE");
        this->mpPhases[4] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[4] == nullptr) << "PhaseTransformationSystem:: UPPER_BAINITE is not set!" << std::endl;
        this->mpPhases[4]->Initialize(rMaterialProperties, rPhaseType);
    }
    if(rMaterialProperties.Has(PHASE_LOWER_BAINITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_LOWER_BAINITE");
        this->mpPhases[5] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[5] == nullptr) << "PhaseTransformationSystem:: LOWER_BAINITE is not set!" << std::endl;
        this->mpPhases[5]->Initialize(rMaterialProperties, rPhaseType);
    }
    if(rMaterialProperties.Has(PHASE_MARTENSITE))
    {
        auto& r_var = KratosComponents<Variable<Phase::Pointer>>::Get("PHASE_MARTENSITE");
        this->mpPhases[6] = rMaterialProperties.GetValue(r_var);
        KRATOS_ERROR_IF(this->mpPhases[6] == nullptr) << "PhaseTransformationSystem:: MARTENSITE is not set!" << std::endl;
        this->mpPhases[6]->Initialize(rMaterialProperties, rPhaseType);
    }

    this->mstat = Status::Set;

    KRATOS_CATCH("");

    iret = 0x1;

    return iret;
}

//************************************************************************************
//************************************************************************************
/**
 * 返回7种相的质量分数
 */
int PhaseTransformationSystem::GetMassFractions(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                              SizeType NodeId, std::array<double, NumPhases>& MassFrac)
{
    int iret = 0x0;
    if(NodeId > 0)
    {
        MassFrac[0] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE);
        MassFrac[1] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_FERRITE);
        MassFrac[2] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_CEMENTITE);
        MassFrac[3] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_PEARLITE);
        MassFrac[4] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE);
        MassFrac[5] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE);
        MassFrac[6] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_MARTENSITE);
        iret = 0x1;
    }
    return iret;
}
int PhaseTransformationSystem::GetLastMassFractions(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                              SizeType NodeId, std::array<double, NumPhases>& LastMassFrac)
{
    int iret = 0x0;
    if(NodeId > 0)
    {
        LastMassFrac[0] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE, 1);
        LastMassFrac[1] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_FERRITE, 1);
        LastMassFrac[2] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_CEMENTITE, 1);
        LastMassFrac[3] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_PEARLITE, 1);
        LastMassFrac[4] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE, 1);
        LastMassFrac[5] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE, 1);
        LastMassFrac[6] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_MARTENSITE, 1);
        iret = 0x1;
    }
    return iret;
}

/**
 * 返回材料的密度
 */
double PhaseTransformationSystem::GetDensity(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0 && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetDensity(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的热导
 */
double PhaseTransformationSystem::GetConductivity(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0 && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetConductivity(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的比热
 */
double PhaseTransformationSystem::GetSpecificHeat(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0 && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetSpecificHeat(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的热导平均密度
 */
double PhaseTransformationSystem::GetDensityMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetDensityMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的平均热导
 */
double PhaseTransformationSystem::GetConductivityMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetConductivityMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的平均比热
 */
double PhaseTransformationSystem::GetSpecificHeatMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetSpecificHeatMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的扩散系数
 */
double PhaseTransformationSystem::GetDC(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetDC(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的平均扩散系数
 */
double PhaseTransformationSystem::GetDCMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetDCMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的杨氏模量
 */
double PhaseTransformationSystem::GetElasticityModulus(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetElasticityModulus(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的泊松比
 */
double PhaseTransformationSystem::GetPoisson(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetPoisson(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的热膨胀系数
 */
double PhaseTransformationSystem::GetExpansion(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetExpansion(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的平均杨氏模量
 */
double PhaseTransformationSystem::GetElasticityModulusMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetElasticityModulus(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的平均泊松比
 */
double PhaseTransformationSystem::GetPoissonMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetPoissonMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的平均热膨胀系数
 */
double PhaseTransformationSystem::GetExpansionMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    GetMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);
    
    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetExpansionMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]);
        }
    }
    return dret;
}

/**
 * 返回材料的屈服应力
 */
double PhaseTransformationSystem::GetYieldStrength(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetConductivity(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的硬化率
 */
double PhaseTransformationSystem::GetHardenRate(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    SizeType i;
    double dret = 0.0;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
            dret += this->mpPhases[i]->GetHardenRate(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
        }
    }
    return dret;
}

/**
 * 返回材料的屈服应力和硬化率
 */
int PhaseTransformationSystem::GetYSandHR(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *ys, double *hr)
{
    int iret;
    SizeType i;
    double dret = 0.0;
    double nys = 0.0, dys;
    double nhr = 0.0, dhr;
    std::array<double, NumPhases> MassFrac;
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);
    
    for(i=0; i<NumPhases; ++i) 
    {
        if(this->mpPhases[i] != nullptr && MassFrac[i] > 0 && MassFrac[i] > 0) {
            iret += this->mpPhases[i]->GetYSandHR(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId, &dys, &dhr);
            nys += dys * MassFrac[i];
            nhr += dhr * MassFrac[i];
        }
    }

    if(iret > 0) 
    {
        if(ys != nullptr) *ys = nys;
        if(hr != nullptr) *hr = nhr;
        iret = 0x1;
    }
    return iret;
}

//************************************************************************************
//************************************************************************************
/**
 * 返回热应变增量
 */
double PhaseTransformationSystem::GetThermalStrainInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                           const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                           SizeType NodeId)
{
    SizeType i;
    double dret = 0.0, d_T = 0.0;
    std::array<double, NumPhases> MassFrac1;
    std::array<double, NumPhases> MassFrac2;

    d_T = rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) - 
              rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1);
    GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac1);
    GetLastMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac2);

    for(i=0; i<NumPhases; ++i) {
        if(MassFrac1[i] > 0 || MassFrac2[i] > 0 && this->mpPhases[i] != nullptr) {
            dret += this->mpPhases[i]->GetHardenRate(rMaterialProperties, rElementGeometry,
                rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * 0.5 * (MassFrac1[i]+MassFrac2[i]) * d_T;
        }
    }
}

/**
 * 返回相变应变增量
 * mfinc[NumPhases]：各相分数增量 
 * mfinc[0]:austenite increment,
 * mfinc[1,2,3] = eutectoid phase incrememnt, ferrtite, cementite, pearlite,
 * mfinc[4,5] = bainite increment, upper bainite, lower bainite,
 * mfinc[6] = martensite increment.
 */
double PhaseTransformationSystem::GetPhaseTransStrainInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *mfinc)
{
    double dret = 0.0;

    if(mfinc != nullptr)
    {
        if(mfinc[0] > 0 && this->mpPhases[0] != nullptr) {
            dret = this->mpPhaseTrans[0]->GetPhaseTransStrain(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * mfinc[0];
        }
        else{
            if(this->mpPhases[1] != nullptr)
            {
                dret = this->mpPhaseTrans[1]->GetPhaseTransStrain(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * (mfinc[1]+mfinc[2]+mfinc[3]);
            }
            if(this->mpPhases[2] != nullptr)
            {
                dret += this->mpPhaseTrans[2]->GetPhaseTransStrain(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * (mfinc[4]+mfinc[5]);
            }
            if(this->mpPhases[3] != nullptr)
            {
                dret += this->mpPhaseTrans[3]->GetPhaseTransStrain(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * mfinc[6];
            }
        }
    }
    return dret;
}

/**
 * 返回蠕变参数
 */
int PhaseTransformationSystem::GetCreepParameters(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                          const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                                          SizeType NodeId, unsigned int *pnm, double *pdata)
{
    int iret = 0x0;
    SizeType i, j, pn;
    double dret = 0.0;
    double dd[10], dd2[10];
    std::array<double, NumPhases> MassFrac;

    
    if(pnm != nullptr && pdata != nullptr)
    {
        GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);

        for(i=0; i<10; ++i) dd[i] = 0.0;

        for(i=0; i<NumPhases; ++i) 
        {
            if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
                if(this->mpPhases[i]->GetCreepParameters(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId, &pn, dd2))
                {
                    for(j=0; j<pn; ++j) dd[j] = dd2[j] * MassFrac[i];
                }
                else
                {
                    break;
                }
            }
        }
        if(i==NumPhases)
        {
            *pnm = pn;
            for(j=0; j<pn; ++i) pdata[j] = dd[j];
            iret = 0x1;
        }
    }
    return iret;
}

/**
 * 返回蠕变应变增量
 */
int PhaseTransformationSystem::GetCreepStrainInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                                        const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                                        SizeType NodeId, double se, double *s, double *cpinc)
{
    int iret = 0x0;
    SizeType i;
    double eqcr = 0.0;
    double dt;

    dt = rCurrentProcessInfo[DELTA_TIME];

    std::array<double, NumPhases> MassFrac;
    
    if(dt > 0 && se >0 && s != nullptr && cpinc != nullptr)
    {
        GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);

        for(i=0; i<NumPhases; ++i) 
        {
            if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
                eqcr += this->mpPhases[i]->GetEqCreepRate(rMaterialProperties, rElementGeometry,
                                                          rShapeFunctionsValues, rCurrentProcessInfo, NodeId) *
                        MassFrac[i];
            }
        }
        eqcr *= 1.5 * dt / se;

        for (i = 0; i < 6; ++i)
        {
            cpinc[i] = eqcr * s[i];
        }

        iret = 0x1;
    }
    
    return iret;
}

/**
 * 返回蠕变应变增量平均值
 */
int PhaseTransformationSystem::GetCreepStrainIncMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double se, double *s, double *cpinc)
{
    int iret = 0x0;
    SizeType i;
    double eqcr = 0.0;
    double dt;

    dt = rCurrentProcessInfo[DELTA_TIME];

    std::array<double, NumPhases> MassFrac;
    
    if(dt > 0 && se >0 && s != nullptr && cpinc != nullptr)
    {
        GetMassFractions(rMaterialProperties, rElementGeometry,
                     rShapeFunctionsValues, rCurrentProcessInfo, NodeId, MassFrac);

        for(i=0; i<NumPhases; ++i) 
        {
            if(this->mpPhases[i] != nullptr && MassFrac[i] > 0) {
                eqcr += this->mpPhases[i]->GetEqCreepRateMean(rMaterialProperties, rElementGeometry,
                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * MassFrac[i];
            }
        }
        eqcr *= 1.5 * dt / se;

        for (i = 0; i < 6; ++i)
        {
            cpinc[i] = eqcr * s[i];
        }

        iret = 0x1;
    }
    
    return iret;
}

/**
 * 返回相变塑性应变增量
 */
int PhaseTransformationSystem::GetTripInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                  const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                                  SizeType NodeId, double *mf, double *mfinc, double *s, double *tpinc)
{
    int iret = 0x0;
	SizeType i;
	double eqtp = 0;
	double af, sf, sfinc;

	if (mf != nullptr && mfinc != nullptr && s != nullptr && tpinc != nullptr)
	{
		if (mfinc[0] > 0 && this->mpPhaseTrans[0] != nullptr)
		{
			eqtp = this->mpPhaseTrans[0]->GetEqTripInc(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId, mf[0], mfinc[0]);
		}
		else
		{
			sf = mf[1] + mf[2] + mf[3];
			sfinc = mfinc[1] + mfinc[2] + mfinc[3];
			af = sf + mf[0];

			if (af > 0 && this->mpPhaseTrans[1] != nullptr)
			{
				sf = sf / af;
				sfinc = sfinc / af;
				eqtp += this->mpPhaseTrans[1]->GetEqTripInc(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId, sf, sfinc) * af;
			}

			sf = mf[4] + mf[5];
			sfinc = mfinc[4] + mfinc[5];
			af = sf + mf[0];

			if (af > 0 && this->mpPhaseTrans[2] != nullptr)
			{
				sf = sf / af;
				sfinc = sfinc / af;
				eqtp += this->mpPhaseTrans[2]->GetEqTripInc(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId, sf, sfinc) * af;
			}

			af = mf[6] + mf[0];
			
			if (af > 0 && this->mpPhaseTrans[3] != nullptr)
			{
				sf = mf[6] / af;
				sfinc = mfinc[6] / af;
				eqtp += this->mpPhaseTrans[3]->GetEqTripInc(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId, sf, sfinc) * af;
			}
		}

		/* get the trip increment */
		for (i = 0; i < 3; i++)
		{
			tpinc[i] = eqtp * s[i];
		}

		for (eqtp *= 2; i < 6; i++)
		{
			tpinc[i] = eqtp * s[i];
		}

		iret = 0x1;
	}

	return iret;
}

/**
 * 返回材料的相变潜热
 */
double PhaseTransformationSystem::GetLatent(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *mfinc)
{
    double dret = 0;

    if (mfinc != nullptr)
    {
        if (mfinc[0] > 0 && this->mpPhaseTrans[0] != nullptr)
        {
            dret = this->mpPhaseTrans[0]->GetPhaseTransLatent(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * mfinc[0];
        }
        else
        {
            if (this->mpPhaseTrans[1] != nullptr)
                dret = this->mpPhaseTrans[1]->GetPhaseTransLatent(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * (mfinc[1] + mfinc[2] + mfinc[3]);
            if (this->mpPhaseTrans[2] != nullptr)
                dret += this->mpPhaseTrans[2]->GetPhaseTransLatent(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * (mfinc[4] + mfinc[5]);
            if (this->mpPhaseTrans[3] != nullptr)
                dret += this->mpPhaseTrans[3]->GetPhaseTransLatent(rMaterialProperties, rElementGeometry,
                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId) * mfinc[6];
        }
    }

    return dret;
}

int PhaseTransformationSystem::GetAusteniteMassFrcInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                              const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                                              SizeType NodeId, double *cfinc, double *mfinc)
{
    int iret = 0x0;
	unsigned int i;
    double t, dt;
	double mt, upt, lwt;
	double mtf;
	double ncf, ncfinc = 0;
	double nmf, nmfinc = 0;
    double stf = 0;

    t = rCurrentProcessInfo[TIME];
    dt = rCurrentProcessInfo[DELTA_TIME];

	if (NodeId > 0 && t >= 0 && dt >= 0 && cfinc != nullptr && mfinc != nullptr)
	{
		*cfinc = 0;
		*mfinc = 0;

		if (this->mpPhaseTrans[0] != nullptr)
		{
			/* get temperature*/
			mt = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                        rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1));
			upt = this->mpPhaseTrans[0]->GetPhaseTransUpperT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
			lwt = this->mpPhaseTrans[0]->GetPhaseTransLowerT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

			if (mt >= lwt && mt <= upt)
			{
				mtf = this->mpPhaseTrans[0]->GetMaxTransMassFrc(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

				if (mtf > 0)
				{
					nmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE, 1);
                    stf = 1.0;

                    ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_AUSTENIZATION, 1);

					/* calculate the austenite mass fraction increment */
					iret = this->mpPhaseTrans[0]->GetPhaseTransMassFrcInc(rMaterialProperties, rElementGeometry,
                                rShapeFunctionsValues, rCurrentProcessInfo, NodeId, ncf, nmf, stf, ncfinc, nmfinc);

					if (iret != 0x0)
					{
						/* adjust the austenite mass fraction increment */
						*mfinc = nmfinc;
						*cfinc = ncfinc;
					}
				}
			}
		}
	}

	/* return the flag */
	return iret;
}

/**
 * 返回共析分解相分数增量
 * mfinc[0]:ferritite
 * mfinc[1]:cementite
 * mfinc[2]:pearlite
 */
int PhaseTransformationSystem::GetEutectoidPhaseMassFrcInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc)
{
    int iret = 0x0;
	unsigned int i;
    double t, dt;
	double mt, upt, lwt;
	double mtf, stf;
	double ncf, ncfinc = 0;
	double nmf, nmfinc = 0;
    double amf, fmf, cmf, pmf;
    double partmf[3];

    t = rCurrentProcessInfo[TIME];
    dt = rCurrentProcessInfo[DELTA_TIME];

    if (NodeId > 0 && t >= 0 && dt >= 0 && cfinc != nullptr && mfinc != nullptr)
    {
        *cfinc = 0;
        for (i = 0; i < 3; i++)     mfinc[i] = 0;

        if (this->mpPhaseTrans[1] != nullptr)
        {
            /* get temperature*/
            mt = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                        rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1));
            upt = this->mpPhaseTrans[1]->GetPhaseTransUpperT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
            lwt = this->mpPhaseTrans[1]->GetPhaseTransLowerT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

            if (mt >= lwt && mt <= upt)
            {
                mtf = this->mpPhaseTrans[1]->GetMaxTransMassFrc(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

                if (mtf > 0)
                {
                    amf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE, 1);
                    fmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_FERRITE, 1);
                    cmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_CEMENTITE, 1);
                    pmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_PEARLITE, 1);
                    stf = amf + fmf + cmf + pmf;
                    nmf = (fmf + cmf + pmf) / stf;

                    ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_EUTECTOID_DECOMPOSITION, 1);

                    /* calculate the eutectoid phase mass fraction increment */
                    iret = this->mpPhaseTrans[1]->GetPhaseTransMassFrcInc(rMaterialProperties, rElementGeometry,
                                rShapeFunctionsValues, rCurrentProcessInfo, NodeId, ncf, nmf, stf, ncfinc, nmfinc);

                    if (iret == 0X1)
                    {
                        /* adjust the eutectoid phase mass fraction increment */
                        /* check mass fraction increment */
						if (nmfinc > 0)
						{
							/* adjust the mass fraction increment */
							nmfinc *= stf;

							/* get partition information */
							iret *= this->mpPhaseTrans[1]->GetPartitionInfo(rMaterialProperties, rElementGeometry,
                                rShapeFunctionsValues, rCurrentProcessInfo, NodeId, amf, partmf);

							/* check partition */
							if (iret == 0x1)
							{
								/* check the partmf[0] */
								if (partmf[0] > 0 && nmfinc > 0)
								{
									if (nmfinc < partmf[0])
									{
										mfinc[0] = nmfinc;
										nmfinc = 0;
									}
									else
									{
										mfinc[0] = partmf[0];
										nmfinc -= partmf[0];
									}
								}
								/* check the partmf[1] */
								if (partmf[1] > 0 && nmfinc > 0)
								{
									if (nmfinc < partmf[1])
									{
										mfinc[1] = nmfinc;
										nmfinc = 0;
									}
									else
									{
										mfinc[1] = partmf[1];
										nmfinc -= partmf[1];
									}
								}
								/* check the partmf[2] */
								if (partmf[2] > 0 && nmfinc > 0)
								{
									if (nmfinc < partmf[2])
									{
										mfinc[2] = nmfinc;
									}
									else
									{
										mfinc[2] = partmf[2];
									}
								}
							}
							else
							{
								/* without partition */
								mfinc[2] = nmfinc;
							}
						}
						/* return the incubation rate */
						*cfinc = ncfinc;
                    }
                }
            }
        }
    }
    return iret;
}

/**
 * 返回贝氏体转变相分数增量
 * mfinc[0]:upper bainite
 * mfinc[1]:lower bainite
 */
int PhaseTransformationSystem::GetBainiteMassFrcInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc)
{
    int iret = 0x0;
	unsigned int i;
    double t, dt;
	double mt, upt, lwt;
	double mtf, stf;
	double ncf, ncfinc = 0.0;
	double nmf, nmfinc = 0.0;
    double amf, ubmf, lbmf;
    double partmf[2];

    t = rCurrentProcessInfo[TIME];
    dt = rCurrentProcessInfo[DELTA_TIME];

    if (NodeId > 0 && t >= 0 && dt >= 0 && cfinc != nullptr && mfinc != nullptr)
    {
        *cfinc = 0.0;
        mfinc[0] = 0.0;
        mfinc[1] = 0.0;

        if (this->mpPhaseTrans[2] != nullptr)
        {
            /* get temperature*/
            mt = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                        rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1));
            upt = this->mpPhaseTrans[2]->GetPhaseTransUpperT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
            lwt = this->mpPhaseTrans[2]->GetPhaseTransLowerT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

            if (mt >= lwt && mt <= upt)
            {
                mtf = this->mpPhaseTrans[2]->GetMaxTransMassFrc(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

                if (mtf > 0)
                {
                    amf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE, 1);
                    ubmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE, 1);
                    lbmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE, 1);
                    stf = amf + ubmf + lbmf;
                    nmf = (ubmf + lbmf) / stf;

                    ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_BAINITIC_TRANSFORMATION, 1);

                    /* calculate the eutectoid phase mass fraction increment */
                    iret = this->mpPhaseTrans[2]->GetPhaseTransMassFrcInc(rMaterialProperties, rElementGeometry,
                                rShapeFunctionsValues, rCurrentProcessInfo, NodeId, ncf, nmf, stf, ncfinc, nmfinc);

                    if (iret == 0X1)
                    {
                        /* adjust the eutectoid phase mass fraction increment */
                        /* check mass fraction increment */
						if (nmfinc > 0)
						{
							/* adjust the mass fraction increment */
							nmfinc *= stf;

							/* get partition information */
							iret *= this->mpPhaseTrans[2]->GetPartitionInfo(rMaterialProperties, rElementGeometry,
                                rShapeFunctionsValues, rCurrentProcessInfo, NodeId, amf, partmf);

							/* check partition */
							if (iret == 0x1)
							{
								/* check the partmf[0] */
								if (partmf[0] > 0 && nmfinc > 0)
								{
									if (nmfinc < partmf[0])
									{
										mfinc[0] = nmfinc;
										nmfinc = 0;
									}
									else
									{
										mfinc[0] = partmf[0];
										nmfinc -= partmf[0];
									}
								}
								/* check the partmf[1] */
								if (partmf[1] > 0 && nmfinc > 0)
								{
									if (nmfinc < partmf[1])
									{
										mfinc[1] = nmfinc;
									}
									else
									{
										mfinc[1] = partmf[1];
									}
								}
							}
							else
							{
								/* without partition */
								mfinc[0] = nmfinc;
							}
						}
						/* return the incubation rate */
						*cfinc = ncfinc;
                    }
                }
            }
        }
    }
    return iret;
}

/**
 * 返回贝氏体转变相分数增量
 */
int PhaseTransformationSystem::GetMartensiteMassFrcInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId, double *cfinc, double *mfinc)
{
    int iret = 0x0;
	unsigned int i;
    double t, dt;
	double mt, upt, lwt;
	double mtf, stf;
	double ncf, ncfinc = 0.0;
	double nmf, nmfinc = 0.0;
    double amf, mmf;

    t = rCurrentProcessInfo[TIME];
    dt = rCurrentProcessInfo[DELTA_TIME];

    if (NodeId > 0 && t >= 0 && dt >= 0 && cfinc != nullptr && mfinc != nullptr)
    {
        *cfinc = 0.0;
        *mfinc = 0.0;

        if (this->mpPhaseTrans[3] != nullptr)
        {
            /* get temperature*/
            mt = rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE);
            upt = this->mpPhaseTrans[3]->GetPhaseTransUpperT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
            lwt = this->mpPhaseTrans[3]->GetPhaseTransLowerT(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

            if (mt >= lwt && mt <= upt)
            {
                mtf = this->mpPhaseTrans[3]->GetMaxTransMassFrc(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);

                if (mtf > 0)
                {
                    amf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE, 1);
                    mmf = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_MARTENSITE, 1);
                    stf = amf + mmf;
                    nmf = mmf / stf;

                    ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_MARTENSITIC_TRANSFORMATION, 1);

                    /* calculate the eutectoid phase mass fraction increment */
                    iret = this->mpPhaseTrans[3]->GetPhaseTransMassFrcInc(rMaterialProperties, rElementGeometry,
                                rShapeFunctionsValues, rCurrentProcessInfo, NodeId, ncf, nmf, stf, ncfinc, nmfinc);

                    if (iret == 0X1)
                    {
                        *mfinc = nmfinc * stf;
                        *cfinc = ncfinc;
                    }
                }
            }
        }
    }
    return iret;
}

/**
 * 更新所有相的质量分数
 */
int PhaseTransformationSystem::UpdateMassFraction(const Properties& rMaterialProperties, GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc)
{
    int iret = 0x0;
	unsigned int i, j;
	double dd, dd2, dd3;
	double ncf, rr;
	double nmf[7];
	double pmf[7];

	/* check the arguments */
	if (NodeId > 0 && cfinc != nullptr && mfinc != nullptr)
	{
		/* part1: update the incubation fraction */
		/* check and update the incubation of austenization */
		if (cfinc[0] > 0)
		{
			ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_AUSTENIZATION, 1) + cfinc[0];
			rElementGeometry[NodeId].GetSolutionStepValue(INCUBATION_FRACTION_AUSTENIZATION) = ncf;
		}
		
		/* check and update the incubation of eutectoid decomposition */
		if (cfinc[1] > 0)
		{
			ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_EUTECTOID_DECOMPOSITION, 1) + cfinc[1];
			rElementGeometry[NodeId].GetSolutionStepValue(INCUBATION_FRACTION_EUTECTOID_DECOMPOSITION) = ncf;
		}

		/* check and update the incubation of bainitic transformation */
		if (cfinc[2] > 0)
		{
			ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_BAINITIC_TRANSFORMATION, 1) + cfinc[1];
			rElementGeometry[NodeId].GetSolutionStepValue(INCUBATION_FRACTION_BAINITIC_TRANSFORMATION) = ncf;
		}

		/* check and update the incubation of martensitic transformation */
		if (cfinc[3] > 0)
		{
			ncf = rElementGeometry[NodeId].FastGetSolutionStepValue(INCUBATION_FRACTION_MARTENSITIC_TRANSFORMATION, 1) + cfinc[1];
			rElementGeometry[NodeId].GetSolutionStepValue(INCUBATION_FRACTION_MARTENSITIC_TRANSFORMATION) = ncf;
		}
		
		/* part2: update the mass fraction */
		/* get the and mass fraction of all microstructures */
		pmf[1] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE, 1);
		pmf[2] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_FERRITE, 1);
		pmf[3] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_CEMENTITE, 1);
		pmf[4] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_PEARLITE, 1);
		pmf[5] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE, 1);
		pmf[6] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE, 1);
		pmf[7] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_MARTENSITE, 1);

		for (i = 0; i < 7; i++)
		{
			nmf[i] = pmf[i];
		}

		if (mfinc[0] > 0)
		{
			/* update the mass fraction of austenite */
			dd = 1 - mfinc[0] / (1 - nmf[0]);
			nmf[0] += mfinc[0];
			nmf[1] *= dd;
			nmf[2] *= dd;
			nmf[3] *= dd;
			nmf[4] *= dd;
			nmf[5] *= dd;
			nmf[6] *= dd;

            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_AUSTENITE) = nmf[0];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_FERRITE) = nmf[1];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_CEMENTITE) = nmf[2];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_PEARLITE) = nmf[3];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE) = nmf[4];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE) = nmf[5];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_MARTENSITE) = nmf[6];

			/* update the carbon content of austenite */
            double cc1, cc2 = 0;

            cc1 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_PEARLITE, 1);
			dd = (pmf[2] - nmf[2]) * 6.69 + (pmf[3] - nmf[3]) * cc1;

            cc1 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_UPPER_BAINITE, 1);
            cc2 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE, 1);
			dd += (pmf[4] - nmf[4]) * cc1 + (pmf[5] - nmf[5]) * cc2;

			cc1 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_MARTENSITE, 1);
			dd += (pmf[6] - nmf[6]) * cc1;

			cc1 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_AUSTENITE, 1);
			dd += pmf[0] * cc1;
			dd /= nmf[0];

            rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_AUSTENITE) = dd;

			/* set the complete flag */
			iret = 0x1;
		}
		else
		{
			/* get the total mass fraction increment */
			for (dd = 0, i = 1; i < 7; i++)
			{
				dd += mfinc[i];
			}

			/* check the total mass fraction increment */
			if (dd > nmf[0])
			{
				rr = nmf[0] / dd;
				dd = nmf[0];
			}
			else
			{
				rr = 1;
			}

			/* update the mass fraction of all microstructures, except austenite */
			for (i = 1; i < 7; i++)
			{
				mfinc[i] *= rr;
				nmf[i] += mfinc[i];
			}

			/* get the new mass fraction of austenite */
			nmf[0] -= dd;

			/* update the mass fraction of all microstructures */
			rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_AUSTENITE) = nmf[0];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_FERRITE) = nmf[1];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_CEMENTITE) = nmf[2];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_PEARLITE) = nmf[3];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE) = nmf[4];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE) = nmf[5];
            rElementGeometry[NodeId].GetSolutionStepValue(MASS_FRACTION_MARTENSITE) = nmf[6];

			/* update the carbon content of austenite */
			dd3 = pmf[0] - mfinc[1] - mfinc[2];

            double cc3;
            cc3 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_AUSTENITE, 1);
			if (dd3 > 0)
			{		
				dd = (pmf[0] * cc3 - mfinc[2] * 6.69) / dd3;
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_AUSTENITE) = dd;
			}
			else
			{
				dd = 0;
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_AUSTENITE) = dd;
			}		

			/* update the carbon content of pearlite */
            cc3 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_PEARLITE, 1);
			if (nmf[3] > 0)
			{
				dd2 = (pmf[3] * cc3 + mfinc[3] * dd) / nmf[3];
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_PEARLITE) = dd2;
			}
			else
			{
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_PEARLITE) = 0;
			}	

			/* update the carbon content of bainite */
            cc3 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_UPPER_BAINITE, 1);
			if (nmf[4] > 0)
			{
				dd2 = (pmf[4] * cc3 + mfinc[4] * dd) / nmf[4];
                rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_UPPER_BAINITE) = dd2;
			}
			else
			{
                rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_UPPER_BAINITE) = 0;
			}

			cc3 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE, 1);
			if (nmf[5] > 0)
			{
				dd2 = (pmf[5] * cc3 + mfinc[5] * dd) / nmf[5];
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE) = dd2;
			}
			else
			{
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE) = 0;
			}

			/* update the carbon content of martensite */
			cc3 = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_MARTENSITE, 1);
			if (nmf[6] > 0)
			{
				dd2 = (pmf[6] * cc3 + mfinc[6] * dd) / nmf[6];
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE) = dd2;
			}
			else
			{
				rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE) = 0;
			}
			/* set the complete flag */
			iret = 0x1;
		}
	}
	/* return the flag */
	return iret;
}

/**
 * 更新所有相的碳含量
 */
int PhaseTransformationSystem::UpdateCarbonContent(const Properties& rMaterialProperties, GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId)
{
    int iret = 0x0;
	unsigned int i;
	double pmf[7], pcc[5];
	double mcc, dd;

	/* check the arguments */
	if (NodeId > 0)
	{
		/* update the carbon content according to the new mean content */
        pmf[1] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_AUSTENITE);
		pmf[2] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_FERRITE);
		pmf[3] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_CEMENTITE);
		pmf[4] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_PEARLITE);
		pmf[5] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_UPPER_BAINITE);
		pmf[6] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_LOWER_BAINITE);
		pmf[7] = rElementGeometry[NodeId].FastGetSolutionStepValue(MASS_FRACTION_MARTENSITE);

        pcc[1] = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_AUSTENITE);
		pcc[4] = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_PEARLITE);
		pcc[5] = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_UPPER_BAINITE);
		pcc[6] = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_LOWER_BAINITE);
		pcc[7] = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_MARTENSITE);

		mcc = rElementGeometry[NodeId].FastGetSolutionStepValue(MEAN_CARBON_CONTENT);
		dd = pmf[2] * 6.69 + pmf[3] * pcc[1] + pmf[4] * pcc[2] + pmf[5] * pcc[3] + pmf[6] * pcc[4];

		/* check the carbon content of cementite and the mass fraction of austenite */		
		if (dd <= mcc && pmf[0] > 0)
		{
			/* set the new carbon content of austenite */
            rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_AUSTENITE) = (mcc - dd) / pmf[0];
		}
		else
		{
			/* rewrite the mean content */
            rElementGeometry[NodeId].GetSolutionStepValue(MEAN_CARBON_CONTENT) = dd;

			/* update the carbon content of austenite */
            rElementGeometry[NodeId].GetSolutionStepValue(CARBON_CONTENT_AUSTENITE) = 0;
		}
		/* set the complete flag */
		iret = 0x1;
	}
	/* return the flag */
	return iret;
}

//************************************************************************************
//************************************************************************************

} // namespace Kratos