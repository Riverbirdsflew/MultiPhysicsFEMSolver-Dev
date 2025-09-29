//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Author:   whf

#include "therm_diff_struct_phase_components.h"
#include "includes/kratos_components.h"

// 你的模块头文件
#include "phase_transformation/kinetics/kinetics_model.h"
#include "phase_transformation/partition/partition_model.h"
#include "phase_transformation/trip/trip_model.h"
#include "phase_transformation/creep/creep_model.h"

namespace Kratos {

// 为相变模块专门定义的AddComponent函数

void AddKratosComponent(const std::string& rName, const KineticsModel& rComponent)
{
    KratosComponents<KineticsModel>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const PartitionModel& rComponent)
{
    KratosComponents<PartitionModel>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const TripModel& rComponent)
{
    KratosComponents<TripModel>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const CreepModel& rComponent)
{
    KratosComponents<CreepModel>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<KineticsModel::Pointer>& rComponent)
{
    KratosComponents<Variable<KineticsModel::Pointer>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<PartitionModel::Pointer>& rComponent)
{
    KratosComponents<Variable<PartitionModel::Pointer>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<TripModel::Pointer>& rComponent)
{
    KratosComponents<Variable<TripModel::Pointer>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<CreepModel::Pointer>& rComponent)
{
    KratosComponents<Variable<CreepModel::Pointer>>::Add(rName, rComponent);
}

// 显式模板实例化
template class KratosComponents<KineticsModel>;
template class KratosComponents<Variable<KineticsModel::Pointer>>;
template class KratosComponents<PartitionModel>;
template class KratosComponents<Variable<PartitionModel::Pointer>>;
template class KratosComponents<TripModel>;
template class KratosComponents<Variable<TripModel::Pointer>>;
template class KratosComponents<CreepModel>;
template class KratosComponents<Variable<CreepModel::Pointer>>;

} // namespace Kratos