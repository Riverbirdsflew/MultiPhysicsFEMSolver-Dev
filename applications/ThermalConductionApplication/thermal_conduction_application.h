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

#if !defined(KRATOS_KRATOS_THERMAL_CONDUCTION_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_THERMAL_CONDUCTION_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "thermal_conduction_application_variables.h"

#include "custom_elements/axisymmetric_eulerian_thermal_conduction.h"
#include "custom_elements/therm_cond_2d.h"
#include "custom_elements/therm_cond_3d.h"
#include "custom_elements/eulerian_cond.h"
#include "custom_elements/eulerian_therm_cond.h"
#include "custom_elements/laplacian_element.h"

#include "custom_conditions/axisymmetric_thermal_face.h"
#include "custom_conditions/thermal_face.h"
#include "custom_conditions/flux_condition.h"

#include "includes/variables.h"
#include "includes/condition.h"


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
class KRATOS_API(THERMAL_CONDUCTION_APPLICATION) KratosThermalConductionApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosThermalConductionApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosThermalConductionApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosThermalConductionApplication();
    
    /// Destructor.
    virtual ~KratosThermalConductionApplication() {}

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

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
        return "KratosThermalConductionApplication";
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

    const AxisymmetricEulerianThermalConductionElement<2,3>  mAxisyEulerianThermCond2D3N;
    const AxisymmetricEulerianThermalConductionElement<2,4>  mAxisyEulerianThermCond2D4N;

    const EulerianThermalConductionElement<2,3>  mEulerianThermCond2D3N;
    const EulerianThermalConductionElement<2,4>  mEulerianThermCond2D4N;
    const EulerianThermalConductionElement<3,4>  mEulerianThermCond3D4N;
    const EulerianThermalConductionElement<3,8>  mEulerianThermCond3D8N;
    const EulerianConductionElement<2,3>  mEulerianCond2D3N;
    const EulerianConductionElement<3,4>  mEulerianCond3D4N;
    const ThermCond2D  mThermCond2D;
    const ThermCond3D  mThermCond3D;
    const LaplacianElement mLaplacian2D3N;
    const LaplacianElement mLaplacian2D4N;
    const LaplacianElement mLaplacian3D4N;
    const LaplacianElement mLaplacian3D8N;

    const AxisymmetricThermalFace mAxisyThermalFace2D2N;
    const ThermalFace mThermalFace2D2N;
    const ThermalFace mThermalFace3D3N;
    const ThermalFace mThermalFace3D4N;
    const FluxCondition<2>  mFluxCondition2D2N;
    const FluxCondition<3>  mFluxCondition3D3N;
    const FluxCondition<4>  mFluxCondition3D4N;

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

    /// Assignment operator.
    KratosThermalConductionApplication& operator=(KratosThermalConductionApplication const& rOther);

    /// Copy constructor.
    KratosThermalConductionApplication(KratosThermalConductionApplication const& rOther);

}; // Class KratosThermalConductionApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_KRATOS_THERMAL_CONDUCTION_APPLICATION_H_INCLUDED  defined
