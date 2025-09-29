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
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "thermal_conduction_application.h"
#include "includes/variables.h"


namespace Kratos {

KratosThermalConductionApplication::KratosThermalConductionApplication():
    KratosApplication("ThermalConductionApplication"),
    
      mAxisyEulerianThermCond2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mAxisyEulerianThermCond2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianThermCond2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianThermCond2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianThermCond3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianThermCond3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))),
      mEulerianCond2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianCond3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mThermCond2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mThermCond3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mLaplacian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))),
      
      mAxisyThermalFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mThermalFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mThermalFace3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mThermalFace3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mFluxCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mFluxCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mFluxCondition3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node >(Element::GeometryType::PointsArrayType(4))))
    {}

void KratosThermalConductionApplication::Register()
{
     KRATOS_INFO("") << "Initializing KratosThermalConductionApplication..." << std::endl;

    //相变潜热
    KRATOS_REGISTER_VARIABLE(PHASE_TRANSFORMATION_LATENT)
    
    // Registering variables
    KRATOS_REGISTER_VARIABLE(AUX_FLUX)
    KRATOS_REGISTER_VARIABLE(AUX_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_1)
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_2)
    KRATOS_REGISTER_VARIABLE(BFECC_ERROR)
    KRATOS_REGISTER_VARIABLE(BFECC_ERROR_1)

    KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
    KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
    KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
    KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)

    KRATOS_REGISTER_VARIABLE(TRANSFER_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(ADJOINT_HEAT_TRANSFER)
    KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
    
    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("AxisyEulerianThermCondElement2D3N", mAxisyEulerianThermCond2D3N);
    KRATOS_REGISTER_ELEMENT("AxisyEulerianThermCondElement2D4N", mAxisyEulerianThermCond2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianThermCondElement2D", mEulerianThermCond2D3N); //TODO: To be removed as it does not follow the naming convention
    KRATOS_REGISTER_ELEMENT("EulerianThermCondElement2D3N", mEulerianThermCond2D3N);
    KRATOS_REGISTER_ELEMENT("EulerianThermCondElement2D4N", mEulerianThermCond2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianThermCondElement3D", mEulerianThermCond3D4N); //TODO: To be removed as it does not follow the naming convention
    KRATOS_REGISTER_ELEMENT("EulerianThermCondElement3D4N", mEulerianThermCond3D4N);
    KRATOS_REGISTER_ELEMENT("EulerianThermCondElement3D8N", mEulerianThermCond3D8N);
    KRATOS_REGISTER_ELEMENT("EulerianCondElement2D3N", mEulerianCond2D3N);
    KRATOS_REGISTER_ELEMENT("EulerianCondElement3D4N", mEulerianCond3D4N);
    KRATOS_REGISTER_ELEMENT("ThermCond2D", mThermCond2D);
    KRATOS_REGISTER_ELEMENT("ThermCond3D", mThermCond3D);
    KRATOS_REGISTER_ELEMENT("LaplacianElement2D3N", mLaplacian2D3N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement2D4N", mLaplacian2D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D4N", mLaplacian3D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D8N", mLaplacian3D8N);

    KRATOS_REGISTER_CONDITION("AxisyThermalFace2D2N", mAxisyThermalFace2D2N);
    KRATOS_REGISTER_CONDITION("ThermalFace2D2N", mThermalFace2D2N);
    KRATOS_REGISTER_CONDITION("ThermalFace3D3N", mThermalFace3D3N);
    KRATOS_REGISTER_CONDITION("ThermalFace3D4N", mThermalFace3D4N);
    KRATOS_REGISTER_CONDITION("FluxCondition2D2N", mFluxCondition2D2N);
    KRATOS_REGISTER_CONDITION("FluxCondition3D3N", mFluxCondition3D3N);
    KRATOS_REGISTER_CONDITION("FluxCondition3D4N", mFluxCondition3D4N);

}

}  // namespace Kratos.
