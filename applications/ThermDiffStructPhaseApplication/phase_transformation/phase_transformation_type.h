//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    WHF
//

#ifndef KRATOS_PHASE_TRANSFORMATION_TYPE_H_INCLUDED
#define KRATOS_PHASE_TRANSFORMATION_TYPE_H_INCLUDED

#include <string>
#include <unordered_map>

namespace Kratos
{
    /**
     * @brief 相变类型枚举
     * 定义不同的相变动力学类型，用于识别使用的动力学模型类型
     */
    enum class PhaseTransformationType
    {
        UNDEFINED_TRANSFORMATION = -1, // 未设定
        AUSTENIZATION = 0,             // 奥氏体转变
        EUTECTOID_DECOMPOSITION = 1,   // 共析转变
        BAINITIC_TRANSFORMATION = 2,    // 贝氏体转变
        MARTENSITIC_TRANSFORMATION = 3, // 马氏体转变
        NUM_TRANSFORMATIONS = 4
    };

    /**
     * @brief 物质相类型枚举
     * 定义钢铁材料中的七种主要相态
     */
    enum class PhaseType
    {
        UNDEFINED = -1,     // 未定义
        AUSTENITE = 0,      // 奥氏体
        FERRITE = 1,        // 铁素体
        CEMENTITE = 2,      // 渗碳体
        PEARLITE = 3,       // 珠光体
        UPPER_BAINITE = 4,  // 上贝氏体
        LOWER_BAINITE = 5,  // 下贝氏体
        MARTENSITE = 6,     // 马氏体
        NUM_PHASES = 7,     //相的种类数
        NO_PHASE = 8       // 没有相成分
    };

    /**
     * @brief 从字符串获取相变类型
     * @param type_name 相变类型的字符串名称
     * @return 对应的PhaseTransformationType枚举值
     */
    inline PhaseTransformationType GetPhaseTransTypeFromString(const std::string& type_name)
    {
        static const std::unordered_map<std::string, PhaseTransformationType> transformation_map = {
            // 英文名称
            {"undefined_transformation", PhaseTransformationType::UNDEFINED_TRANSFORMATION},
            {"austenization", PhaseTransformationType::AUSTENIZATION},
            {"eutectoid_decomposition", PhaseTransformationType::EUTECTOID_DECOMPOSITION},
            {"bainitic_transformation", PhaseTransformationType::BAINITIC_TRANSFORMATION},
            {"martensitic_transformation", PhaseTransformationType::MARTENSITIC_TRANSFORMATION},
            
            // 简化名称
            {"undefined", PhaseTransformationType::UNDEFINED_TRANSFORMATION},
            {"austenite", PhaseTransformationType::AUSTENIZATION},
            {"eutectoid", PhaseTransformationType::EUTECTOID_DECOMPOSITION},
            {"bainitic", PhaseTransformationType::BAINITIC_TRANSFORMATION},
            {"bainite", PhaseTransformationType::BAINITIC_TRANSFORMATION},
            {"martensitic", PhaseTransformationType::MARTENSITIC_TRANSFORMATION},
            {"martensite", PhaseTransformationType::MARTENSITIC_TRANSFORMATION},
            
            // 数字字符串
            {"-1", PhaseTransformationType::UNDEFINED_TRANSFORMATION},
            {"0", PhaseTransformationType::AUSTENIZATION},
            {"1", PhaseTransformationType::EUTECTOID_DECOMPOSITION},
            {"2", PhaseTransformationType::BAINITIC_TRANSFORMATION},
            {"3", PhaseTransformationType::MARTENSITIC_TRANSFORMATION}
        };

        auto it = transformation_map.find(type_name);
        return (it != transformation_map.end()) ? it->second : PhaseTransformationType::UNDEFINED_TRANSFORMATION;
    }

    /**
     * @brief 从字符串获取相态类型
     * @param type_name 相态类型的字符串名称
     * @return 对应的PhaseType枚举值
     */
    inline PhaseType GetPhaseTypeFromString(const std::string& type_name)
    {
        static const std::unordered_map<std::string, PhaseType> phase_map = {
            // 英文名称
            {"undefined", PhaseType::UNDEFINED},
            {"austenite", PhaseType::AUSTENITE},
            {"ferrite", PhaseType::FERRITE},
            {"cementite", PhaseType::CEMENTITE},
            {"pearlite", PhaseType::PEARLITE},
            {"upper_bainite", PhaseType::UPPER_BAINITE},
            {"lower_bainite", PhaseType::LOWER_BAINITE},
            {"martensite", PhaseType::MARTENSITE},
            {"no_phase", PhaseType::NO_PHASE},
            
            // 数字字符串
            {"-1", PhaseType::UNDEFINED},
            {"0", PhaseType::AUSTENITE},
            {"1", PhaseType::FERRITE},
            {"2", PhaseType::CEMENTITE},
            {"3", PhaseType::PEARLITE},
            {"4", PhaseType::UPPER_BAINITE},
            {"5", PhaseType::LOWER_BAINITE},
            {"6", PhaseType::MARTENSITE},
            {"8", PhaseType::NO_PHASE}
        };

        auto it = phase_map.find(type_name);
        return (it != phase_map.end()) ? it->second : PhaseType::UNDEFINED;
    }

} // namespace Kratos

#endif // KRATOS_PHASE_DEFINITIONS_H_INCLUDED