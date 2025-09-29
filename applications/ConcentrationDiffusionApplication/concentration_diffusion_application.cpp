//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoffi @{KRATOS_APP_AUTHOR}
//


// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/line_2d_2.h"
#include "geometries/point_2d.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "concentration_diffusion_application.h"
#include "concentration_diffusion_application_variables.h"
#include "includes/variables.h"


namespace Kratos {

KratosConcentrationDiffusionApplication::KratosConcentrationDiffusionApplication():
    KratosApplication("ConcentrationDiffusionApplication"),
      //轴对称单元
      mAxisymmetricEulerianConcentrationDiffusion2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymmetricEulerianConcentrationDiffusion2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      //三角形、四边形、四面体、六面体单元
      mEulerianConcenDiff2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianConcenDiff2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConcenDiff3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConcenDiff3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))),
      mEulerianDiffusion2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianDiffusion3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mConcenDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mConcenDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      
      //轴对称面对流边界
      mAxisymmetricConcentrationFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      //平面、三维面对流边界
      mConcentrationFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mConcentrationFace3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mConcentrationFace3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      //第二类边界
      mConcentrationFluxCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mConcentrationFluxCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mConcentrationFluxCondition3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node >(Element::GeometryType::PointsArrayType(4))))
    {}

void KratosConcentrationDiffusionApplication::Register()
{
     //KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosConcentrationDiffusionApplication..." << std::endl;
  
  // Registering variables//注册变量
  KRATOS_REGISTER_VARIABLE( AMBIENT_CONCENTRATION )
  KRATOS_REGISTER_VARIABLE( FACE_CONCENTRATION_FLUX )
  KRATOS_REGISTER_VARIABLE( CONCENTRATION_FLUX )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONCENTRATION_GRADIENT)

  KRATOS_REGISTER_VARIABLE( PHASE_TRANSITION_FRACTION )
  KRATOS_REGISTER_VARIABLE( PHASE_NAME )
  
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

  KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
  KRATOS_REGISTER_VARIABLE(TRANSFER_COEFFICIENT)
  
    
    //~ KRATOS_REGISTER_VARIABLE( CONCENTRATION_FACE )
    //~ KRATOS_REGISTER_VARIABLE(AUX_FLUX)
    //~ KRATOS_REGISTER_VARIABLE(AUX_TEMPERATURE)
    //~ KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_1)
    //~ KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_2)
    //~ KRATOS_REGISTER_VARIABLE(BFECC_ERROR)
    //~ KRATOS_REGISTER_VARIABLE(BFECC_ERROR_1)

    //~ KRATOS_REGISTER_VARIABLE(MEAN_SIZE)

    //~ KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
    //~ KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)


    //~ KRATOS_REGISTER_VARIABLE(ADJOINT_HEAT_TRANSFER)
    //~ KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTION)


  
  // Registering elements and conditions here//注册单元
    KRATOS_REGISTER_ELEMENT("AxisyEulerianConcenDiff2D3NElement", mAxisymmetricEulerianConcentrationDiffusion2D3N);
    KRATOS_REGISTER_ELEMENT("AxisyEulerianConcenDiff2D4NElement", mAxisymmetricEulerianConcentrationDiffusion2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConcenDiff2D", mEulerianConcenDiff2D3N); //TODO: To be removed as it does not follow the naming convention
    KRATOS_REGISTER_ELEMENT("EulerianConcenDiff2D3NElement", mEulerianConcenDiff2D3N);
    KRATOS_REGISTER_ELEMENT("EulerianConcenDiff2D4NElement", mEulerianConcenDiff2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConcenDiff3D", mEulerianConcenDiff3D4N); //TODO: To be removed as it does not follow the naming convention
    KRATOS_REGISTER_ELEMENT("EulerianConcenDiff3D4NElement", mEulerianConcenDiff3D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConcenDiff3D8NElement", mEulerianConcenDiff3D8N);
    KRATOS_REGISTER_ELEMENT("EulerianDiffusion2D3NElement", mEulerianDiffusion2D3N);
    KRATOS_REGISTER_ELEMENT("EulerianDiffusion3D4NElement", mEulerianDiffusion3D4N);
    KRATOS_REGISTER_ELEMENT("ConcenDiff2D", mConcenDiff2D);
    KRATOS_REGISTER_ELEMENT("ConcenDiff3D", mConcenDiff3D);
   
    //轴对称面对流边界
    KRATOS_REGISTER_CONDITION("AxisymmetricConcentrationFace2D2N", mAxisymmetricConcentrationFace2D2N);
    //面对流边界
    KRATOS_REGISTER_CONDITION("ConcentrationFace2D2N", mConcentrationFace2D2N);
    KRATOS_REGISTER_CONDITION("ConcentrationFace3D3N", mConcentrationFace3D3N);
    KRATOS_REGISTER_CONDITION("ConcentrationFace3D4N", mConcentrationFace3D4N);
    //第二类边界
    KRATOS_REGISTER_CONDITION("ConcentrationFluxCondition2D2N", mConcentrationFluxCondition2D2N);
    KRATOS_REGISTER_CONDITION("ConcentrationFluxCondition3D3N", mConcentrationFluxCondition3D3N);
    KRATOS_REGISTER_CONDITION("ConcentrationFluxCondition3D4N", mConcentrationFluxCondition3D4N);

}

}  // namespace Kratos.
