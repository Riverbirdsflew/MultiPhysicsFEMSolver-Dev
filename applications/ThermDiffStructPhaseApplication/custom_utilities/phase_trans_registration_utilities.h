// includes/registration_utilities.h
#ifndef PHASE_TRANS_REGISTRATION_UTILITIES_H_INCLUDED
#define PHASE_TRANS_REGISTRATION_UTILITIES_H_INCLUDED

// Kratos Core includes
#include "includes/kratos_components.h"
#include "includes/registry.h"
#include "includes/serializer.h"

// 模块相关的include
#include "phase_transformation/creep/creep_model.h"
#include "phase_transformation/kinetics/kinetics_model.h"
#include "phase_transformation/partition/partition_model.h"
#include "phase_transformation/trip/trip_model.h"

namespace Kratos {

/**
 * @brief 通用组件注册模板函数（接受智能指针）
 * @param name 组件名称
 * @param reference 组件对象指针
 * @param component_type_name 组件类型名称（用于Registry键值）
 */
template<typename TComponentType>
void RegisterComponent(const std::string& name, 
                      const TComponentType& reference,
                      const std::string& component_type_name) 
{
    // 添加到KratosComponents
    KratosComponents<TComponentType>::Add(name, reference);
    
    // 构造Registry键值
    std::string registry_key = component_type_name + "." + Registry::GetCurrentSource() + "." + name;
    std::string components_key = "components." + name;
    
    // 注册到Registry（避免重复注册）
    if (!Registry::HasItem(registry_key) && !Registry::HasItem(components_key)) {
        Registry::AddItem<RegistryItem>(registry_key);
        Registry::AddItem<RegistryItem>(components_key);
    }
    
    // 序列化注册
    Serializer::Register(name, reference);
}

// /**
//  * @brief 通用组件注册模板函数（接受对象引用，会创建shared_ptr）
//  * @param name 组件名称
//  * @param reference 组件对象引用
//  * @param component_type_name 组件类型名称（用于Registry键值）
//  */
// template<typename TComponentType>
// void RegisterComponent(const std::string& name, 
//                       const TComponentType& reference,
//                       const std::string& component_type_name) 
// {
//     // 创建shared_ptr包装对象
//     auto ptr_reference = Kratos::make_shared<TComponentType>(reference);
    
//     // 调用指针版本的注册函数
//     RegisterComponent<TComponentType>(name, ptr_reference, component_type_name);
// }

/**
 * @brief CreepModel专用注册函数
 * @param name 模型名称
 * @param reference CreepModel对象指针
 */
inline void RegisterCreepModel(const std::string& name, const CreepModel& reference) {
    RegisterComponent<CreepModel>(name, reference, "creep_model");
}

/**
 * @brief KineticsModel专用注册函数
 * @param name 模型名称  
 * @param reference KineticsModel对象指针
 */
inline void RegisterKineticsModel(const std::string& name, const KineticsModel& reference) {
    RegisterComponent<KineticsModel>(name, reference, "kinetics_model");
}

/**
 * @brief PartitionModel专用注册函数
 * @param name 模型名称
 * @param reference PartitionModel对象指针
 */
inline void RegisterPartitionModel(const std::string& name, const PartitionModel& reference) {
    RegisterComponent<PartitionModel>(name, reference, "partition_model");
}

/**
 * @brief TripModel专用注册函数
 * @param name 模型名称
 * @param reference TripModel对象指针  
 */
inline void RegisterTripModel(const std::string& name, const TripModel& reference) {
    RegisterComponent<TripModel>(name, reference, "trip_model");
}

} // namespace Kratos

#endif // PHASE_TRANS_REGISTRATION_UTILITIES_H_INCLUDED