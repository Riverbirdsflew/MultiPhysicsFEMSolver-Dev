//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#pragma once

// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"
#include "phase_transformation/creep/creep_model.h"
#include "phase_transformation/kinetics/kinetics_model.h"
#include "phase_transformation/trip/trip_model.h"
#include "phase_transformation/partition/partition_model.h"
#include "phase_transformation/phase.h"
#include "phase_transformation/phase_transformation_law.h"


namespace Kratos {

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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) KratosThermDiffStructPhaseApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosThermDiffStructPhaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosThermDiffStructPhaseApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosThermDiffStructPhaseApplication();

    /// Copy constructor.
    KratosThermDiffStructPhaseApplication(KratosThermDiffStructPhaseApplication const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    KratosThermDiffStructPhaseApplication& operator=(KratosThermDiffStructPhaseApplication const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    /// 重写注销方法来处理相变的组件
    void DeregisterApplication() override;

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
    std::string Info() const override
    {
        return "KratosThermDiffStructPhaseApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
          KRATOS_WATCH("in my application");
          KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    // Base phase transformation definition
    const CreepModel mCreepModel;
    const KineticsModel mKineticsModel;
    const PartitionModel mPartitionModel;
    const TripModel mTripModel;
    const Phase mPhase;
    const PhaseTransformationLaw mPhaseTransformationLaw;


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

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    // const Elem2D   mElem2D;
    // const Elem3D   mElem3D;

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

}; // Class KratosThermDiffStructPhaseApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
