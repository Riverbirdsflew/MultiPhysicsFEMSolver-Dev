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

# pragma once

#ifndef THERM_DIFF_STRUCT_PHASE_COMPONENTS_H_INCLUDED
#define THERM_DIFF_STRUCT_PHASE_COMPONENTS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos {

// 前向声明
class KineticsModel;
class PartitionModel; 
class TripModel;
class CreepModel;
template<class TDataType> class Variable;


// 声明相变模块的组件添加函数
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const KineticsModel& rComponent);
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const PartitionModel& rComponent);
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const TripModel& rComponent);
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const CreepModel& rComponent);

void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const Variable<KineticsModel::Pointer>& rComponent);
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const Variable<PartitionModel::Pointer>& rComponent);
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const Variable<TripModel::Pointer>& rComponent);
void KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) AddKratosComponent(const std::string& rName, const Variable<CreepModel::Pointer>& rComponent);

} // namespace Kratos

#endif // THERM_DIFF_STRUCT_PHASE_COMPONENTS_H_INCLUDED