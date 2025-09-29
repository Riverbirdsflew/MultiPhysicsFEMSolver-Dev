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

// System includes

// External includes

// Project includes
#include "utilities/read_materials_utility.h"
#include "includes/kratos_components.h"

namespace Kratos {

/**
 * @class ReadPhaseTransMaterialsUtility
 * @ingroup KratosCore
 * @brief Read phase-transformation-related models from json and assign to Properties
 *
 * JSON 支持四个相变过程：
 *   - austenitization
 *   - eutectoid
 *   - bainitic
 *   - martensitic
 * 对每个过程可配置：kinetics / partition / trip 三类模型。
 *
 * 同时支持七种相的蠕变模型：
 *   - austenite, ferrite, cementite, pearlite, martensite, upper_bainite, lower_bainite
 *
 * 变量（需在 Application 中注册）：
 *   - Variable<KineticsModel::Pointer>:
 *       AUSTENITIZATION_KINETICS_MODEL, EUTECTOID_KINETICS_MODEL,
 *       BAINITIC_KINETICS_MODEL, MARTENSITIC_KINETICS_MODEL
 *   - Variable<PartitionModel::Pointer>:
 *       AUSTENITIZATION_PARTITION_MODEL, EUTECTOID_PARTITION_MODEL,
 *       BAINITIC_PARTITION_MODEL, MARTENSITIC_PARTITION_MODEL
 *   - Variable<TripModel::Pointer>:
 *       AUSTENITIZATION_TRIP_MODEL, EUTECTOID_TRIP_MODEL,
 *       BAINITIC_TRIP_MODEL, MARTENSITIC_TRIP_MODEL
 *   - Variable<CreepModel::Pointer>:
 *       AUSTENITE_CREEP_MODEL, FERRITE_CREEP_MODEL, CEMENTITE_CREEP_MODEL,
 *       PEARLITE_CREEP_MODEL, MARTENSITE_CREEP_MODEL,
 *       UPPER_BAINITE_CREEP_MODEL, LOWER_BAINITE_CREEP_MODEL
 */

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

///@}
///@name Kratos Classes
///@{

/**
 * @class ReadMaterialsUtility
 * @ingroup KratosCore
 * @brief Process to read constitutive law and material properties from a json file
 * @details This process reads constitutive law and material properties from a json file
 * and assign them to elements and conditions.
 * The definition includes the creation of subproperties
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) ReadPhaseTransMaterialsUtility : public ReadMaterialsUtility
{
  public:

    ///@name Type Definitions
    ///@{

    using BaseType = ReadMaterialsUtility;

      ///@}
      ///@name Pointer Definitions

      /// Pointer definition of ReadMaterialProcess
      KRATOS_CLASS_POINTER_DEFINITION(ReadPhaseTransMaterialsUtility);

      ///@}
      ///@name Life Cycle
      ///@{
      
      /**
       * @brief Default constructor
       * @param rModel The model containing the problem to solve
       */
      ReadPhaseTransMaterialsUtility(Model &rModel) : BaseType(rModel) {};

      /**
       * @brief Constructor reading directly from file, via parameters
       * @param Params The configuration parameters telling where the configuration file can be found
       * @param rModel The model containing the problem to solve
       */
      ReadPhaseTransMaterialsUtility(Parameters Params,Model &rModel) : BaseType(Params, rModel) {}

      /**
       * @brief Constructor reading directly from file, via text
       * @param Params The string telling where the configuration file can be found
       * @param rModel The model containing the problem to solve
       */
      ReadPhaseTransMaterialsUtility(const std::string &rParametersName, Model &rModel) : BaseType(rParametersName, rModel) {}

      /// Destructor.
      virtual ~ReadPhaseTransMaterialsUtility() override = default;

      ///@}
      ///@name Operators
      ///@{

      ///@}
      ///@name Operations
      ///@{

      /**
       * @brief This reads the properties from parameters
       * @param MaterialData The configuration parameters defining the properties
       */
      void ReadMaterials(Parameters MaterialData);

      // === Overrides ===
      /** 扩展材料赋值：先使用基类，再解析相变/蠕变配置 */
      void AssignMaterialToProperty(
          const Parameters MaterialData,
          Properties &rProperty) override;

      /**
       解析四类相变模型（动力学/配分/TRIP）
       */
      void AssignPhaseTransLawToProperty(
          const Parameters MaterialData,
          Properties &rProperty);

        /**
        // 解析七种相的蠕变模型
       */
      void AssignCreepToProperty(
          const Parameters MaterialData,
          Properties &rProperty);


      ///@}
      ///@name Access
      ///@{

      ///@}
      ///@name Inquiry
      ///@{

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      std::string Info() const
      {
          return "ReadPhaseTransMaterialsUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const  {
        rOStream << "ReadPhaseTransMaterialsUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const  {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class ReadMaterialsUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
