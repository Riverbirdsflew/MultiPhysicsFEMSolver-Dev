// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| ANALYSIS
//
//  License:         BSD License
//                   license: StructuralAnalysisApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @namespace RayleighDampingCoefficientsUtilities
 * @ingroup StructuralAnalysisApplication
 * @brief This utility computes the two first eigen values of the system and estimates the alpha and beta Rayleigh damping coefficients // 计算系统的前两个特征值,估计雷利阻尼系数,有效地描述结构在动态加载下的能量耗散特性，在地震和振动分析中
 * @details It uses the well stablished formulation: Wilson, E. L. (2004). Static and Dynamic Analysis of Structures (4th ed.). Berkeley, CA: Computers and Structures, Inc.
 * @author Vicente Mataix Ferrandiz
*/
namespace RayleighDampingCoefficientsUtilities
{
    /**
    * @brief This utility computes the two first eigen values of the system and estimates the alpha and beta Rayleigh damping coefficients
    * @details It uses the well stablished formulation: Wilson, E. L. (2004). Static and Dynamic Analysis of Structures (4th ed.). Berkeley, CA: Computers and Structures, Inc.
    * @param ThisParameters The configuration parameters
    */
    Vector KRATOS_API(STRUCTURAL_ANALYSIS_APPLICATION) ComputeDampingCoefficients(Parameters ThisParameters);

} /// namespace RayleighDampingCoefficientsUtilities
} /// namespace Kratos
