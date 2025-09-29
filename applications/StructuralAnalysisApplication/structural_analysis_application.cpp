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


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "structural_analysis_application.h"
#include "structural_analysis_application_variables.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/pyramid_3d_5.h"
#include "geometries/pyramid_3d_13.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

namespace Kratos {
	
	// define node type
	typedef Node NodeType;

KratosStructuralAnalysisApplication::KratosStructuralAnalysisApplication():
    KratosApplication("StructuralAnalysisApplication"),
    /* ELEMENTS */
    // Adding the nodal concentrated element
    mNodalConcentratedElement2D1N(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))), true),
    mNodalConcentratedDampedElement2D1N(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))), false),
    mNodalConcentratedElement3D1N(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))), true),
    mNodalConcentratedDampedElement3D1N(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))), false),    

    // Adding the kinematic linear elements
      mSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mSmallDisplacement2D10N(0, Element::GeometryType::Pointer(new Triangle2D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mSmallDisplacement2D15N(0, Element::GeometryType::Pointer(new Triangle2D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mSmallDisplacement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacement3D5N(0, Element::GeometryType::Pointer(new Pyramid3D5<NodeType >(Element::GeometryType::PointsArrayType(5)))),
      mSmallDisplacement3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacement3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mSmallDisplacement3D13N(0, Element::GeometryType::Pointer(new Pyramid3D13<NodeType >(Element::GeometryType::PointsArrayType(13)))),
      mSmallDisplacement3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mSmallDisplacement3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mSmallDisplacement3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),
      
      mAxisymSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAxisymSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAxisymSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),

      mSmallDisplacementMixedVolumetricStrainElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementMixedVolumetricStrainElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementMixedVolumetricStrainElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementMixedVolumetricStrainElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      
      /* CONDITIONS */
      // Adding point load conditions
      mPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAxisymPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      // Adding line load conditions
      mLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mLineLoadCondition2D4N(0, Condition::GeometryType::Pointer(new Line2D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mLineLoadCondition2D5N(0, Condition::GeometryType::Pointer(new Line2D5<NodeType >(Condition::GeometryType::PointsArrayType(5)))),
      mLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mLineLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Line3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mSmallDisplacementLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementLineLoadCondition2D4N(0, Condition::GeometryType::Pointer(new Line2D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementLineLoadCondition2D5N(0, Condition::GeometryType::Pointer(new Line2D5<NodeType >(Condition::GeometryType::PointsArrayType(5)))),
      mSmallDisplacementLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mSmallDisplacementLineLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Line3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAxisymLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAxisymLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      // Adding surface load conditions
      mSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >( Condition::GeometryType::PointsArrayType(4)))),
      mSurfaceLoadCondition3D6N(0, Condition::GeometryType::Pointer(new Triangle3D6<NodeType >(Condition::GeometryType::PointsArrayType(6)))),
      mSurfaceLoadCondition3D8N(0, Condition::GeometryType::Pointer(new Quadrilateral3D8<NodeType >(Condition::GeometryType::PointsArrayType(8)))),
      mSurfaceLoadCondition3D9N(0, Condition::GeometryType::Pointer(new Quadrilateral3D9<NodeType >(Condition::GeometryType::PointsArrayType(9)))),
      mSmallDisplacementSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >( Condition::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementSurfaceLoadCondition3D6N(0, Condition::GeometryType::Pointer(new Triangle3D6<NodeType >(Condition::GeometryType::PointsArrayType(6)))),
      mSmallDisplacementSurfaceLoadCondition3D8N(0, Condition::GeometryType::Pointer(new Quadrilateral3D8<NodeType >(Condition::GeometryType::PointsArrayType(8)))),
      mSmallDisplacementSurfaceLoadCondition3D9N(0, Condition::GeometryType::Pointer(new Quadrilateral3D9<NodeType >(Condition::GeometryType::PointsArrayType(9)))),
      
      // Adding the displacement-control condition
      mDisplacementControlCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1))))
    {}

void KratosStructuralAnalysisApplication::Register()
{
    KRATOS_INFO("") << "    KRATOS   ___|  |                   |                   |\n"
                    << "           \\___ \\  __|  __| |   |  __| __| |   |  __| _` | |\n"
                    << "                 | |   |    |   | (    |   |   | |   (   | |\n"
                    << "           _____/ \\__|_|   \\__,_|\\___|\\__|\\__,_|_|  \\__,_|_| ANALYSIS\n"
                    << "Initializing KratosStructuralAnalysisApplication..." << std::endl;

    // General pourpose
    KRATOS_REGISTER_VARIABLE(INTEGRATION_ORDER); // The integration order considered on the element
    KRATOS_REGISTER_VARIABLE(LOCAL_MATERIAL_AXIS_1);
    KRATOS_REGISTER_VARIABLE(LOCAL_MATERIAL_AXIS_2);
    KRATOS_REGISTER_VARIABLE(LOCAL_MATERIAL_AXIS_3);
    KRATOS_REGISTER_VARIABLE(LOCAL_PRESTRESS_AXIS_1);
    KRATOS_REGISTER_VARIABLE(LOCAL_PRESTRESS_AXIS_2);
    KRATOS_REGISTER_VARIABLE(CENTER_OF_GRAVITY);
    KRATOS_REGISTER_VARIABLE(MASS_MOMENT_OF_INERTIA);
    KRATOS_REGISTER_VARIABLE(ELASTICITY_TENSOR);


    // Generalized eigenvalue problem
    KRATOS_REGISTER_VARIABLE(BUILD_LEVEL)
    KRATOS_REGISTER_VARIABLE(EIGENVALUE_VECTOR)
    KRATOS_REGISTER_VARIABLE(EIGENVECTOR_MATRIX)
    KRATOS_REGISTER_VARIABLE(MODAL_MASS_MATRIX)
    KRATOS_REGISTER_VARIABLE(MODAL_STIFFNESS_MATRIX)

    // Energy
    KRATOS_REGISTER_VARIABLE(ENERGY_DAMPING_DISSIPATION)

    // Mixed formulations generalized variables
    KRATOS_REGISTER_VARIABLE( REACTION_STRAIN )

    // Formfinding
    KRATOS_REGISTER_VARIABLE(LAMBDA_MAX)
    KRATOS_REGISTER_VARIABLE(IS_FORMFINDING)
    KRATOS_REGISTER_VARIABLE(BASE_REF_1)
    KRATOS_REGISTER_VARIABLE(BASE_REF_2)

    // Cross section
    KRATOS_REGISTER_VARIABLE(SHELL_CROSS_SECTION_OUTPUT_PLY_ID)
    KRATOS_REGISTER_VARIABLE(SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_LAYERS)

    // Nodal stiffness
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_INITIAL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_DISPLACEMENT_STIFFNESS)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_INITIAL_ROTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_STIFFNESS)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_DAMPING_RATIO)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_DAMPING_RATIO)

    // CONDITIONS
    // Reset equations ids "flag"
    KRATOS_REGISTER_VARIABLE(RESET_EQUATION_IDS);

    // Strain measures
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_VECTOR);
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_TENSOR);

    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS)

    KRATOS_REGISTER_VARIABLE(REFERENCE_DEFORMATION_GRADIENT);
    KRATOS_REGISTER_VARIABLE(REFERENCE_DEFORMATION_GRADIENT_DETERMINANT);

    // Rayleigh variables
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    // System damping
    KRATOS_REGISTER_VARIABLE(SYSTEM_DAMPING_RATIO)
    KRATOS_REGISTER_VARIABLE(SECOND_SYSTEM_DAMPING_RATIO)

    // Nodal load variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MOVING_LOAD)
    KRATOS_REGISTER_VARIABLE(MOVING_LOAD_LOCAL_DISTANCE)

    // Condition load variables
    KRATOS_REGISTER_VARIABLE(POINT_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(LINE_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(MOVING_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(POSITIVE_FACE_PRESSURES_VECTOR)
    KRATOS_REGISTER_VARIABLE(NEGATIVE_FACE_PRESSURES_VECTOR)

    // Displacement-Control variables
    KRATOS_REGISTER_VARIABLE(LOAD_FACTOR)
    KRATOS_REGISTER_VARIABLE(PRESCRIBED_DISPLACEMENT)

    // Register the nodal concentrated element
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement2D1N", mNodalConcentratedElement2D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement3D1N", mNodalConcentratedElement3D1N);

    // SOLID ELEMENTS
    // Small displacement elements
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D3N", mSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D4N", mSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D6N", mSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D8N", mSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D9N", mSmallDisplacement2D9N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D10N", mSmallDisplacement2D10N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D15N", mSmallDisplacement2D15N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D4N", mSmallDisplacement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D5N", mSmallDisplacement3D5N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D6N", mSmallDisplacement3D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D8N", mSmallDisplacement3D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D10N", mSmallDisplacement3D10N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D13N", mSmallDisplacement3D13N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D15N", mSmallDisplacement3D15N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D20N", mSmallDisplacement3D20N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D27N", mSmallDisplacement3D27N)

    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement2D3N", mSmallDisplacementMixedVolumetricStrainElement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement2D4N", mSmallDisplacementMixedVolumetricStrainElement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement3D4N", mSmallDisplacementMixedVolumetricStrainElement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement3D8N", mSmallDisplacementMixedVolumetricStrainElement3D8N)

    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D3N", mAxisymSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D4N", mAxisymSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D6N", mAxisymSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D8N", mAxisymSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D9N", mAxisymSmallDisplacement2D9N)


    // Register the conditions
    // Point loads
    KRATOS_REGISTER_CONDITION("PointLoadCondition2D1N", mPointLoadCondition2D1N)
    KRATOS_REGISTER_CONDITION("PointLoadCondition3D1N", mPointLoadCondition3D1N)

    KRATOS_REGISTER_CONDITION("AxisymPointLoadCondition2D1N", mAxisymPointLoadCondition2D1N)

    // Line loads
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D2N", mLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D3N", mLineLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D4N", mLineLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D5N", mLineLoadCondition2D5N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition3D2N", mLineLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition3D3N", mLineLoadCondition3D3N)

    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D2N", mSmallDisplacementLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D3N", mSmallDisplacementLineLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D4N", mSmallDisplacementLineLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D5N", mSmallDisplacementLineLoadCondition2D5N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition3D2N", mSmallDisplacementLineLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition3D3N", mSmallDisplacementLineLoadCondition3D3N)

    KRATOS_REGISTER_CONDITION("AxisymLineLoadCondition2D2N", mAxisymLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("AxisymLineLoadCondition2D3N", mAxisymLineLoadCondition2D3N)

    // Surface loads
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D6N", mSurfaceLoadCondition3D6N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D8N", mSurfaceLoadCondition3D8N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D9N", mSurfaceLoadCondition3D9N)

    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D3N", mSmallDisplacementSurfaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D4N", mSmallDisplacementSurfaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D6N", mSmallDisplacementSurfaceLoadCondition3D6N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D8N", mSmallDisplacementSurfaceLoadCondition3D8N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D9N", mSmallDisplacementSurfaceLoadCondition3D9N)

    // Displacement-Control Conditions
    KRATOS_REGISTER_CONDITION("DisplacementControlCondition3D1N", mDisplacementControlCondition3D1N)


    // Register linear elastics laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic3DLaw", mElasticIsotropic3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrain2DLaw", mLinearPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStress2DLaw", mLinearPlaneStress);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticAxisym2DLaw", mAxisymElasticIsotropic);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("UserProvidedLinearElastic2DLaw", mUserProvidedLinearElastic2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("UserProvidedLinearElastic3DLaw", mUserProvidedLinearElastic3DLaw);
}

}  // namespace Kratos.
