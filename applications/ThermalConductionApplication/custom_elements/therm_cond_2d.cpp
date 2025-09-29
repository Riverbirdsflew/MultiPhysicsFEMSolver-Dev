// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes 

// External includes 

// Project includes 
#include "includes/define.h"
#include "custom_elements/therm_cond_2d.h"
#include "thermal_conduction_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

ThermCond2D::ThermCond2D(IndexType NewId, GeometryType::Pointer pGeometry)
: Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ThermCond2D::ThermCond2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
: Element(NewId, pGeometry, pProperties)
{

}

//************************************************************************************
//************************************************************************************

Element::Pointer ThermCond2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermCond2D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer ThermCond2D::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermCond2D>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

ThermCond2D::~ThermCond2D()
{
}
  
  //************************************************************************************
  //************************************************************************************
  
  void ThermCond2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
                                        VectorType& rRightHandSideVector, 
                                        const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int number_of_points = GetGeometry().size();
    const double lumping_factor = 1.00 / double(number_of_points);
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    double dt_inv = 1.0 / delta_t;
    unsigned int TDim = 2;
    const int Stationary = rCurrentProcessInfo[STATIONARY];
    if (rLeftHandSideMatrix.size1() != number_of_points)
      rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
    
    if (rRightHandSideVector.size() != number_of_points)
      rRightHandSideVector.resize(number_of_points, false);
    
    BoundedMatrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
    BoundedMatrix<double, 3, 2 > msDN_DX;
    array_1d<double, 3 > msN;
    array_1d<double, 3 > ms_temp_vec_np;
    array_1d<double, 3 > ms_u_DN;
    
    
    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    
    
    const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
    const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
    
    double conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
    //        double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
    double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(rSpecificHeatVar);
    double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
    double heat_flux = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
    
    for (unsigned int i = 1; i < number_of_points; i++)
      {
	conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
	density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
	
	specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
	heat_flux += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);
	
      }
    conductivity *= lumping_factor;
    density *= lumping_factor;
    specific_heat *= lumping_factor;
    heat_flux *= lumping_factor;
    
    //INERTIA CONTRIBUTION
    if(Stationary!=1){
      noalias(rLeftHandSideMatrix) += dt_inv * (density * specific_heat) * msMassFactors;
    }
    // RHS = Fext (heat per unit mass)
    noalias(rRightHandSideVector) = (heat_flux * density) * msN;
    
    //adding the inertia terms
    // RHS += M*vhistory
    //calculating the historical velocity
    
    if(Stationary!=1){
      for (unsigned int iii = 0; iii < number_of_points; iii++)
	ms_temp_vec_np[iii] = -dt_inv * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);
      noalias(rRightHandSideVector) -= prod(msMassFactors, ms_temp_vec_np * density * specific_heat);
    }
    //subtracting the dirichlet term
    // RHS -= LHS*temperatures
    for (unsigned int iii = 0; iii < number_of_points; iii++)
      ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_temp_vec_np);
    
    
    rRightHandSideVector *= Area;
    rLeftHandSideMatrix *= Area;
    
    KRATOS_CATCH("");
  }
  
  

  //************************************************************************************
  //************************************************************************************

    void ThermCond2D::CalculateRightHandSide(VectorType& rRightHandSideVector, 
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void ThermCond2D::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
    {
    }


    //************************************************************************************
    //************************************************************************************

    void ThermCond2D::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const 
    {
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    }

    //************************************************************************************
    //************************************************************************************

    void ThermCond2D::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

    }
 
    //************************************************************************************
    //************************************************************************************

} // Namespace Kratos


