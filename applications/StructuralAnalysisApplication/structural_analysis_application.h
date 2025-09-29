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
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

/* ELEMENTS */
/* 0D elements */
#include "custom_elements/nodal_concentrated_element.hpp"
/* Adding solid elements */
#include "custom_elements/small_displacement.h"
#include "custom_elements/axisym_small_displacement.h"

/* Adding the mixed solid elements */
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"

/* Conditions */
#include "custom_conditions/base_load_condition.h"
#include "custom_conditions/point_load_condition.h"
#include "custom_conditions/axisym_point_load_condition.h"
#include "custom_conditions/line_load_condition.h"
#include "custom_conditions/small_displacement_line_load_condition.h"
#include "custom_conditions/axisym_line_load_condition_2d.h"
#include "custom_conditions/surface_load_condition_3d.h"
#include "custom_conditions/small_displacement_surface_load_condition_3d.h"
#include "custom_conditions/displacement_control_condition.h"

/* Constitutive Laws */
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/user_provided_linear_elastic_law.h"



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

/// Short class definition.KratosStructuralAnalysisApplication
/** Detail class definition.Application for Structural analysis
*/
class KRATOS_API(STRUCTURAL_ANALYSIS_APPLICATION) KratosStructuralAnalysisApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosStructuralAnalysisApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralAnalysisApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosStructuralAnalysisApplication();
    
    ///Deconstructor.
    ~KratosStructuralAnalysisApplication() override {}


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
        return "KratosStructuralAnalysisApplication";
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
        KRATOS_INFO("KratosStructuralAnalysisApplication") << "Has the following number of components: " << KratosComponents<VariableData>::GetComponents().size() << std::endl;

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

    /* ELEMENTS */
    // Adding the nodal concentrated element
    const NodalConcentratedElement mNodalConcentratedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedElement3D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement3D1N;
    // Linear kinematic elements
    const SmallDisplacement mSmallDisplacement2D3N;
    const SmallDisplacement mSmallDisplacement2D4N;
    const SmallDisplacement mSmallDisplacement2D6N;
    const SmallDisplacement mSmallDisplacement2D8N;
    const SmallDisplacement mSmallDisplacement2D9N;
    const SmallDisplacement mSmallDisplacement2D10N;
    const SmallDisplacement mSmallDisplacement2D15N;
    const SmallDisplacement mSmallDisplacement3D4N;
    const SmallDisplacement mSmallDisplacement3D5N;
    const SmallDisplacement mSmallDisplacement3D6N;
    const SmallDisplacement mSmallDisplacement3D8N;
    const SmallDisplacement mSmallDisplacement3D10N;
    const SmallDisplacement mSmallDisplacement3D13N;
    const SmallDisplacement mSmallDisplacement3D15N;
    const SmallDisplacement mSmallDisplacement3D20N;
    const SmallDisplacement mSmallDisplacement3D27N;
     
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D3N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D4N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D6N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D8N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D9N;

    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement2D3N;
    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement2D4N;
    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement3D4N;
    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement3D8N;
    

    /* CONDITIONS*/
    // Point load
    const PointLoadCondition mPointLoadCondition2D1N;
    const PointLoadCondition mPointLoadCondition3D1N;

    const AxisymPointLoadCondition mAxisymPointLoadCondition2D1N;

    // Line load
    const LineLoadCondition<2> mLineLoadCondition2D2N;
    const LineLoadCondition<2> mLineLoadCondition2D3N;
    const LineLoadCondition<2> mLineLoadCondition2D4N;
    const LineLoadCondition<2> mLineLoadCondition2D5N;
    const LineLoadCondition<3> mLineLoadCondition3D2N;
    const LineLoadCondition<3> mLineLoadCondition3D3N;

    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D2N;
    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D3N;
    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D4N;
    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D5N;
    const SmallDisplacementLineLoadCondition<3> mSmallDisplacementLineLoadCondition3D2N;
    const SmallDisplacementLineLoadCondition<3> mSmallDisplacementLineLoadCondition3D3N;

    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D2N;
    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D3N;

    // Surface load
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D3N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D4N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D6N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D8N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D9N;

    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D3N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D4N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D6N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D8N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D9N;
    
    // Displacement-Control Conditions
    const DisplacementControlCondition mDisplacementControlCondition3D1N;
    
    /* CONSTITUTIVE LAWS */
    // Linear elastics laws
    const ElasticIsotropic3D mElasticIsotropic3D;
    const AxisymElasticIsotropic mAxisymElasticIsotropic;
    const LinearPlaneStrain  mLinearPlaneStrain;
    const LinearPlaneStress  mLinearPlaneStress;
    const UserProvidedLinearElasticLaw<2> mUserProvidedLinearElastic2DLaw;
    const UserProvidedLinearElasticLaw<3> mUserProvidedLinearElastic3DLaw;
    
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
    KratosStructuralAnalysisApplication& operator=(KratosStructuralAnalysisApplication const& rOther);

    /// Copy constructor.
    KratosStructuralAnalysisApplication(KratosStructuralAnalysisApplication const& rOther);
    ///@}
    
}; // Class KratosStructuralAnalysisApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
