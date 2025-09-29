// KRATOS ___ ___  _  _  ___ ___ _  _     ___ ___ ___ ___ 
//       / __/ _ \| \| |/ __| __| \| |___|   \_ _| __| __|
//      | (_| (_) | .` | (__| _|| .` |___| |) | || _|| _| 
//       \___\___/|_|\_|\___|___|_|\_|   |___/___|_| |_|  APPLICATION
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
#include "custom_elements/concen_diff_3d.h"
#include "concentration_diffusion_application.h"
#include "includes/concentration_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ConcenDiff3D::ConcenDiff3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ConcenDiff3D::ConcenDiff3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer ConcenDiff3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcenDiff3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer ConcenDiff3D::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcenDiff3D>(NewId, pGeom, pProperties);
}

ConcenDiff3D::~ConcenDiff3D()
{
}

//************************************************************************************
//************************************************************************************

void ConcenDiff3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

      const unsigned int number_of_points = GetGeometry().size();
    const double lumping_factor = 1.00 / double(number_of_points);
    unsigned int TDim = 3;
    const int Stationary = rCurrentProcessInfo[STATIONARY];		
    
    const BoundedMatrix<double, 4, 4 > msMassFactors = 0.25 * IdentityMatrix(4, 4);
    BoundedMatrix<double, 4, 3 > msDN_DX;
    array_1d<double, 4 > msN;
    array_1d<double, 3 > ms_vel_gauss;
    array_1d<double, 4 > ms_concen_vec_np;
    array_1d<double, 4 > ms_u_DN;
    
    BoundedMatrix<double, 3, 3 > First = ZeroMatrix(3, 3);
    BoundedMatrix<double, 3, 3 > Second = ZeroMatrix(3, 3);
    BoundedMatrix<double, 3, 4 > Third = ZeroMatrix(3, 4);
    BoundedMatrix<double, 3, 3 > Identity = 1.0 * IdentityMatrix(3, 3);
    
    array_1d<double, 3 > grad_g = ZeroVector(3); //dimesion coincides with space dimension

    
    if (rLeftHandSideMatrix.size1() != number_of_points)
      rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
    
    if (rRightHandSideVector.size() != number_of_points)
      rRightHandSideVector.resize(number_of_points, false);
    
    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    ConcentrationDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);
    
    //const Variable<double>& rDensityVar = my_settings->GetDensityVariable();//密度
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();
    const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
    //const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();//比热
    
    
    double diffuse_conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
    //double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);//密度
    //double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(rSpecificHeatVar);//比热	
    //double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
    double concen_flux = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);//concen flux物质流
    double proj = GetGeometry()[0].FastGetSolutionStepValue(rProjectionVariable);
    //double nu = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    
    //    const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double, 3 > & v= GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar); 
    const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
    for (unsigned int j = 0; j < TDim; j++)
      ms_vel_gauss[j] = v[j] - w[j];
    
    for (unsigned int i = 1; i < number_of_points; i++)
      {
        diffuse_conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
        //density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);求密度总和
	//specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);//求比热总和
        //specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
        concen_flux += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);//求concen flux物质流总和
        proj += GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable);
        //nu += GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
	
	const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
	//        const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
        for (unsigned int j = 0; j < TDim; j++)
	  ms_vel_gauss[j] += v[j] - w[j];
	
    }
    diffuse_conductivity *= lumping_factor;
    //density *= lumping_factor;//密度平均
    //specific_heat *= lumping_factor;//比热平均
    concen_flux *= lumping_factor;//concen flux物质流密度平均
    proj *= lumping_factor;
    ms_vel_gauss *= lumping_factor;
    //    nu *= lumping_factor;
    
    //getting the BDF2 coefficients (not fixed to allow variable time step)
    //the coefficients INCLUDE the time step
    const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    
    for (unsigned int i = 0; i < TDim; i++)
      {
        for (unsigned int j = 0; j < number_of_points; j++)
	  {
            grad_g[i] += msDN_DX(j, i) * GetGeometry()[j].FastGetSolutionStepValue(rUnknownVar);
	  }
      }
    
    //double norm_g = norm_2(grad_g);

    double res = (inner_prod(ms_vel_gauss, grad_g)); //+ 0.333333333333333 * (t0media+t1media+t2media)*(1/dt)*density*diffuse_conductivity;
    double aux_res = fabs(res - proj);
    if (fabs(res) > aux_res)
      res = aux_res;
    //        res -= proj;
    //res *= density*specific_heat;//密度*比热=1
    double norm_grad = norm_2(grad_g);
    double k_aux = fabs(res) / (norm_grad + 1e-6);
    k_aux *= 0.707;
    
    
    //calculating parameter tau
    double c1 = 4.00;
    double c2 = 2.00;
    double h = pow(6.00 * Area, 0.3333333);
    double norm_u = norm_2(ms_vel_gauss);
    double tau1;
    if(Stationary==1){
      //tau1 = (h * h) / ( c1 * diffuse_conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);//密度*比热=1
      tau1 = (h * h) / ( c1 * diffuse_conductivity + c2 * (norm_u + 1e-6) * h);
    }
    else{
      //tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * diffuse_conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);//密度*比热=1
      tau1 = (h * h) / (BDFcoeffs[0] * h * h + c1 * diffuse_conductivity + c2 *(norm_u + 1e-6) * h);
    }

    noalias(First) = outer_prod(ms_vel_gauss, trans(ms_vel_gauss));
    First /= ((norm_u + 1e-6)*(norm_u + 1e-6));
    noalias(Second) = Identity - First;
    noalias(Third) = prod(Second, trans(msDN_DX));
    
    //CONCENTRATION CONTRIBUTION TO THE STIFFNESS MATRIX
    noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
    //noalias(rLeftHandSideMatrix) = (density * specific_heat) * outer_prod(msN, ms_u_DN);//密度*比热=1
    noalias(rLeftHandSideMatrix) = outer_prod(msN, ms_u_DN);
    
    //CONCENTRATION STABILIZING CONTRIBUTION (Suu)
    
    //noalias(rLeftHandSideMatrix) += density * specific_heat * density * specific_heat * tau1 * outer_prod(ms_u_DN, ms_u_DN);//密度*比热=1
    noalias(rLeftHandSideMatrix) += tau1 * outer_prod(ms_u_DN, ms_u_DN);
    
    //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
    noalias(rLeftHandSideMatrix) += diffuse_conductivity  * prod(msDN_DX, trans(msDN_DX)) + k_aux * h * prod(msDN_DX, Third);
    
    if(Stationary!=1){
      //INERTIA CONTRIBUTION
      //noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * (density * specific_heat) * msMassFactors;//密度*比热=1
      noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMassFactors;
    }
    
    // RHS = Fext
    //noalias(rRightHandSideVector) = (heat_flux * density) * msN;//内热源
    noalias(rRightHandSideVector) = (concen_flux) * msN;
    
    //RHS += Suy * proj[component]
    //noalias(rRightHandSideVector) += density * specific_heat * density * specific_heat * (tau1 * proj) * ms_u_DN;//密度*比热=1
    noalias(rRightHandSideVector) += (tau1 * proj) * ms_u_DN;
    if(Stationary!=1){
      //adding the inertia terms
      // RHS += M*vhistory
      //calculating the historical velocity
      for (unsigned int iii = 0; iii < number_of_points; iii++)
        ms_concen_vec_np[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);
      for (unsigned int step = 2; step < BDFcoeffs.size(); step++)
	{
	  for (unsigned int iii = 0; iii < number_of_points; iii++)
            ms_concen_vec_np[iii] += BDFcoeffs[step] * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, step);
	}
      //noalias(rRightHandSideVector) -= prod(msMassFactors, ms_concen_vec_np * density * specific_heat);//密度*比热=1
      noalias(rRightHandSideVector) -= prod(msMassFactors, ms_concen_vec_np);
    }
    //subtracting the dirichlet term
    // RHS -= LHS*concentration
    for (unsigned int iii = 0; iii < number_of_points; iii++)
      ms_concen_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_concen_vec_np);
    
    //multiplying by area
    rRightHandSideVector *= Area;
    rLeftHandSideMatrix *= Area;
    
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ConcenDiff3D::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ConcenDiff3D::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    BoundedMatrix<double, 4, 3 > msDN_DX;
    array_1d<double, 4 > msN;
    array_1d<double, 3 > ms_vel_gauss;
    array_1d<double, 4 > ms_concen_vec_np;
    array_1d<double, 4 > ms_u_DN;

    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    ConcentrationDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
    const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();

    if (FractionalStepNumber == 2) //calculation of concentration convective projection
    {
        const unsigned int number_of_points = GetGeometry().size();
        const double lumping_factor = 1.00 / double(number_of_points);
        unsigned int TDim = 3;

        //calculating viscosity
        ms_concen_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
        //const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar);
        const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
        for (unsigned int j = 0; j < TDim; j++)
            ms_vel_gauss[j] = v[j] - w[j];

        for (unsigned int i = 1; i < number_of_points; i++)
        {
            ms_concen_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
	    const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
	    //const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
            for (unsigned int j = 0; j < TDim; j++)
                ms_vel_gauss[j] += v[j] - w[j];

        }
        ms_vel_gauss *= lumping_factor;

        //calculating convective auxiliary vector
        noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
        double concen_conv = inner_prod(ms_u_DN, ms_concen_vec_np);
        concen_conv *= Area;

        for (unsigned int i = 0; i < number_of_points; i++)
        {
            GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
            GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable) += lumping_factor*concen_conv;
            ;
        }
    }
    KRATOS_CATCH("");
}


//************************************************************************************
//************************************************************************************

void ConcenDiff3D::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const 
{

    ConcentrationDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
}

//************************************************************************************
//************************************************************************************

void ConcenDiff3D::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const 
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConcentrationDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if (ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

}


} // Namespace Kratos


