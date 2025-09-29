//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hofii @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "concentration_diffusion_application_variables.h"

#include "custom_elements/concen_diff_2d.h"
#include "custom_elements/concen_diff_3d.h"
#include "custom_elements/eulerian_concen_diff.h"
#include "custom_elements/eulerian_diff.h"
#include "custom_elements/axisymmetric_eulerian_concentration_diffusion.h"

#include "custom_conditions/concentration_flux_condition.h"
#include "custom_conditions/concentration_face_condition.h"
#include "custom_conditions/axisymmetric_concentration_face.h"

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
class KRATOS_API(CONCENTRATION_DIFFUSION_APPLICATION) KratosConcentrationDiffusionApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosConcentrationDiffusionApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosConcentrationDiffusionApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConcentrationDiffusionApplication();

    /// Destructor.
    virtual ~KratosConcentrationDiffusionApplication() {}


    ///@}
    ///@name Operators
    ///@{


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
        return "KratosConcentrationDiffusionApplication";
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
        KRATOS_WATCH("in KratosConcentrationDiffusionApplication");
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

    ///@}
    ///@name Member Variables
    ///@{
    
    const AxisymmetricEulerianConcentrationDiffusionElement<2,3>  mAxisymmetricEulerianConcentrationDiffusion2D3N;
    const AxisymmetricEulerianConcentrationDiffusionElement<2,4>  mAxisymmetricEulerianConcentrationDiffusion2D4N;

    const EulerianConcentrationDiffusionElement<2,3>  mEulerianConcenDiff2D3N;
    const EulerianConcentrationDiffusionElement<2,4>  mEulerianConcenDiff2D4N;
    const EulerianConcentrationDiffusionElement<3,4>  mEulerianConcenDiff3D4N;
    const EulerianConcentrationDiffusionElement<3,8>  mEulerianConcenDiff3D8N;
    const EulerianDiffusionElement<2,3>  mEulerianDiffusion2D3N;
    const EulerianDiffusionElement<3,4>  mEulerianDiffusion3D4N;

    const ConcenDiff2D  mConcenDiff2D;
    const ConcenDiff3D  mConcenDiff3D;

    const AxisymmetricConcentrationFace mAxisymmetricConcentrationFace2D2N;

    const ConcentrationFace mConcentrationFace2D2N;
    const ConcentrationFace mConcentrationFace3D3N;
    const ConcentrationFace mConcentrationFace3D4N;
    const ConcentrationFluxCondition<2>  mConcentrationFluxCondition2D2N;
    const ConcentrationFluxCondition<3>  mConcentrationFluxCondition3D3N;
    const ConcentrationFluxCondition<4>  mConcentrationFluxCondition3D4N;


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
    KratosConcentrationDiffusionApplication& operator=(KratosConcentrationDiffusionApplication const& rOther);

    /// Copy constructor.
    KratosConcentrationDiffusionApplication(KratosConcentrationDiffusionApplication const& rOther);


    ///@}

}; // Class KratosConcentrationDiffusionApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_CONCENTRATION_DIFFUSION_APPLICATION_H_INCLUDED  defined
