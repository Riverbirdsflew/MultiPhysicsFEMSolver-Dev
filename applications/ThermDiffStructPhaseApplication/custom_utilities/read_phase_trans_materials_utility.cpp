//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//



// System includes
#include <filesystem>

// External includes

// Project includes
#include "read_phase_trans_materials_utility.h"
#include "includes/kratos_components.h"
#include "../phase_transformation/kinetics/kinetics_model.h"
#include "../phase_transformation/trip/trip_model.h"
#include "../phase_transformation/creep/creep_model.h"
#include "../phase_transformation/partition/partition_model.h"
#include "../phase_transformation/phase_transformation_law.h"
#include "../phase_transformation/phase.h"

namespace Kratos {

void ReadPhaseTransMaterialsUtility::AssignMaterialToProperty(
    const Parameters MaterialData,
    Properties& rProperty)
{
    // 调用基类方法分配普通材料属性
    BaseType::AssignMaterialToProperty(MaterialData, rProperty);

    // 调用自己的方法处理相变
    AssignPhaseTransLawToProperty(MaterialData, rProperty);

    // 调用自己的方法处理蠕变
    AssignCreepToProperty(MaterialData, rProperty);
}

// ----------------------- 解析四类相变模型 -----------------------
void ReadPhaseTransMaterialsUtility::AssignPhaseTransLawToProperty(
    const Parameters MaterialData,
    Properties& rProperty)
{
    if (!MaterialData.Has("phase_transformation"))
        return;

    Parameters phase_trans_para = MaterialData["phase_transformation"];

    auto assign_one_block = [this, &rProperty](
        const std::string& block_name,
        const char* kinetics_var_name,
        const char* partition_var_name,
        const char* trip_var_name,
        PhaseTransformationType transformation_type,
        const char* law_var_name,
        Parameters& parent)
    {
        if (!parent.Has(block_name)) 
            return;

        Parameters blk = parent[block_name];

        KineticsModel::Pointer p_kinetics = nullptr;
        PartitionModel::Pointer p_partition = nullptr;
        TripModel::Pointer p_trip = nullptr;

        // ---- kinetics ----
        if (blk.Has("kinetics")) {
            Parameters kn = blk["kinetics"];
            std::string kname = kn.GetString();
            this->TrimComponentName(kname);

            ///KRATOS_ERROR_IF_NOT(KratosComponents<KineticsModel>::Has(kname))<< "Missing kinetics model: " << kname << std::endl;

            auto& r_var = KratosComponents<Variable<KineticsModel::Pointer>>::Get(kinetics_var_name);
            p_kinetics = KratosComponents<KineticsModel>::Get(kname).Create(kn, rProperty);
            rProperty.SetValue(r_var, p_kinetics);
        }
        else
        {
            auto& r_var = KratosComponents<Variable<KineticsModel::Pointer>>::Get(kinetics_var_name);
            rProperty.SetValue(r_var, nullptr);
        }

        // ---- partition ----
        if (blk.Has("partition")) {
            Parameters pn = blk["partition"];
            std::string pname = pn.GetString();
            this->TrimComponentName(pname);

            // KRATOS_ERROR_IF_NOT(KratosComponents<PartitionModel>::Has(pname))<< "Missing partition model: " << pname << std::endl;

            auto& r_var = KratosComponents<Variable<PartitionModel::Pointer>>::Get(partition_var_name);
            p_partition = KratosComponents<PartitionModel>::Get(pname).Create(pn, rProperty);
            rProperty.SetValue(r_var, p_partition);
        }
        else
        {
            auto& r_var = KratosComponents<Variable<PartitionModel::Pointer>>::Get(partition_var_name);
            rProperty.SetValue(r_var, nullptr);
        }

        // ---- trip ----
        if (blk.Has("trip")) {
            Parameters tn = blk["trip"];
            std::string tname = tn.GetString();
            this->TrimComponentName(tname);

            // KRATOS_ERROR_IF_NOT(KratosComponents<TripModel>::Has(tname))<< "Missing TRIP model: " << tname << std::endl;

            auto& r_var = KratosComponents<Variable<TripModel::Pointer>>::Get(trip_var_name);
            p_trip = KratosComponents<TripModel>::Get(tname).Create(tn, rProperty);
            rProperty.SetValue(r_var, p_trip);
        }
        else
        {
            auto& r_var = KratosComponents<Variable<TripModel::Pointer>>::Get(trip_var_name);
            rProperty.SetValue(r_var, nullptr);
        }

        // 创建 PhaseTransformationLaw 对象，使用带参数的构造函数
        auto p_phase_law = Kratos::make_shared<PhaseTransformationLaw>(
            transformation_type,
            p_kinetics,
            p_partition, 
            p_trip);
        
        // 存储 PhaseTransformationLaw 对象
        auto& phase_law_var = KratosComponents<Variable<PhaseTransformationLaw::Pointer>>::Get(law_var_name);
        rProperty.SetValue(phase_law_var, p_phase_law);
    };

    // 四种相变块
    assign_one_block(
        "austenitization",
        "KINETICS_AUSTENIZATION",
        "",
        "TRIP_AUSTENIZATION",
        PhaseTransformationType::AUSTENIZATION,
        "PHASE_TRANSFORMATION_LAW_AUSTENIZATION",
        phase_trans_para);

    assign_one_block(
        "eutectoid",
        "KINETICS_EUTECTOID_DECOMPOSITION",
        "PARTITION_EUTECTOID_DECOMPOSITION",
        "TRIP_EUTECTOID_DECOMPOSITION",
        PhaseTransformationType::EUTECTOID_DECOMPOSITION,
        "PHASE_TRANSFORMATION_LAW_EUTECTOID_DECOMPOSITION",
        phase_trans_para);

    assign_one_block(
        "bainitic",
        "KINETICS_BAINITIC_TRANSFORMATION",
        "PARTITION_BAINITIC_TRANSFORMATION",
        "TRIP_EUTECTOID_DECOMPOSITION",
        PhaseTransformationType::BAINITIC_TRANSFORMATION,
        "PHASE_TRANSFORMATION_LAW_BAINITIC_TRANSFORMATION",
        phase_trans_para);

    assign_one_block(
        "martensitic",
        "KINETICS_MARTENSITIC_TRANSFORMATION",
        "",
        "TRIP_MARTENSITIC_TRANSFORMATION",
        PhaseTransformationType::MARTENSITIC_TRANSFORMATION,
        "PHASE_TRANSFORMATION_LAW_MARTENSITIC_TRANSFORMATION",
        phase_trans_para);
}

// ----------------------- 解析七种相的蠕变模型 -----------------------
void ReadPhaseTransMaterialsUtility::AssignCreepToProperty(
    const Parameters MaterialData,
    Properties& rProperty)
{
    if (!MaterialData.Has("creep_models"))
        return;

    Parameters creep = MaterialData["creep_models"];

    auto assign_creep = [this, &rProperty, &creep](
                            const char *json_key,
                            const char *creep_var_name,
                            PhaseType phase_type,
                            const char *phase_var_name)
    {
        if (creep.Has(json_key))
        {
            std::string creep_model_name = creep[json_key].GetString();
            this->TrimComponentName(creep_model_name);

            Parameters cm;
            cm.AddString("name", creep_model_name);

            // KRATOS_ERROR_IF_NOT(KratosComponents<CreepModel>::Has(cname))<< "Missing creep model: " << cname << " for phase " << json_key << std::endl;

            auto &r_var = KratosComponents<Variable<CreepModel::Pointer>>::Get(creep_var_name);
            auto p_model = KratosComponents<CreepModel>::Get(creep_model_name).Create(cm, rProperty);
            rProperty.SetValue(r_var, p_model);

            auto p_phase = Kratos::make_shared<Phase>(phase_type, p_model);
            auto &phase_var = KratosComponents<Variable<Phase::Pointer>>::Get(phase_var_name);
            rProperty.SetValue(phase_var, p_phase);
        }
        else
        {
            auto &r_var = KratosComponents<Variable<CreepModel::Pointer>>::Get(creep_var_name);
            rProperty.SetValue(r_var, nullptr);

            auto p_phase = Kratos::make_shared<Phase>(phase_type, nullptr);
            auto &phase_var = KratosComponents<Variable<Phase::Pointer>>::Get(phase_var_name);
            rProperty.SetValue(phase_var, p_phase);
        }
    };

    assign_creep("austenite",     "CREEP_AUSTENITE",     PhaseType::AUSTENITE,     "PHASE_AUSTENITE" );
    assign_creep("ferrite",       "CREEP_FERRITE",       PhaseType::FERRITE,       "PHASE_FERRITE" );
    assign_creep("cementite",     "CREEP_CEMENTITE",     PhaseType::CEMENTITE,     "PHASE_CEMENTITE" );
    assign_creep("pearlite",      "CREEP_PEARLITE",      PhaseType::PEARLITE,      "PHASE_PEARLITE" );
    assign_creep("upper_bainite", "CREEP_UPPER_BAINITE", PhaseType::UPPER_BAINITE, "PHASE_UPPER_BAINITE" );
    assign_creep("lower_bainite", "CREEP_LOWER_BAINITE", PhaseType::LOWER_BAINITE, "PHASE_LOWER_BAINITE" );
    assign_creep("martensite",    "CREEP_MARTENSITE",    PhaseType::MARTENSITE,    "PHASE_MARTENSITE" );
    
}


}  // namespace Kratos.
