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
#include "custom_elements/concen_diff_2d.h"
#include "concentration_diffusion_application.h"
#include "includes/concentration_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

ConcenDiff2D::ConcenDiff2D(IndexType NewId, GeometryType::Pointer pGeometry)
: Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ConcenDiff2D::ConcenDiff2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
: Element(NewId, pGeometry, pProperties)
{

}

//************************************************************************************
//************************************************************************************

Element::Pointer ConcenDiff2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcenDiff2D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer ConcenDiff2D::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcenDiff2D>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

ConcenDiff2D::~ConcenDiff2D()
{
}
  
  //************************************************************************************
  //************************************************************************************
  
  void ConcenDiff2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
                                        VectorType& rRightHandSideVector, 
                                        const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
      //getting the BDF2 coefficients (not fixed to allow variable time step)
      //the coefficients INCLUDE the time step
      const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    
    const unsigned int number_of_points = GetGeometry().size();
    const double lumping_factor = 1.00 / double(number_of_points);//所有节点平均
    unsigned int TDim = 2;
    const int Stationary = rCurrentProcessInfo[STATIONARY];
    if (rLeftHandSideMatrix.size1() != number_of_points)
      rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
    
    if (rRightHandSideVector.size() != number_of_points)
      rRightHandSideVector.resize(number_of_points, false);
    
    BoundedMatrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);//1/3对角阵
    BoundedMatrix<double, 3, 2 > msDN_DX;//形函数导数阵
    array_1d<double, 3 > msN;
    array_1d<double, 2 > ms_vel_gauss;
    array_1d<double, 3 > ms_concen_vec_np;///未知量/浓度向量
    array_1d<double, 3 > ms_u_DN;
    array_1d<double, 2 > grad_g;
    BoundedMatrix<double, 2, 2 > Identity = IdentityMatrix(2, 2);
    BoundedMatrix<double, 2, 2 > First;
    BoundedMatrix<double, 2, 2 > Second;
    BoundedMatrix<double, 2, 3 > Third;
    
    
    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);//形函数导数，形函数值，面积
    
    ConcentrationDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);//
    
    
    //const Variable<double>& rDensityVar = my_settings->GetDensityVariable(); //扩散求解没有密度变量，删除或者密度*比热=1
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();//扩散系数
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();//未知量
    const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();
    const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
    //const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();//比热变量
    
    double diffuse_conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);//扩散系数
    //double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
    //double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(rSpecificHeatVar);//比热
    //double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);//密度
    double concen_flux = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);//边界物质浓度流
    double proj = GetGeometry()[0].FastGetSolutionStepValue(rProjectionVariable);
    //const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY); 
    const array_1d<double, 3 > & v= GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar); 
    const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); 
    
    for (unsigned int j = 0; j < TDim; j++)
      ms_vel_gauss[j] = v[j] - w[j];//速度量
    
    for (unsigned int i = 1; i < number_of_points; i++) //求单元的平均量
      {
	diffuse_conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
	//density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);//所有点密度的加和
	
	//specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);//所有点比热的加和
	//specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
	concen_flux += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);////所有点物质浓度流的加和
	proj += GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable);
	//const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);	
	const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
	for (unsigned int j = 0; j < TDim; j++)
	  ms_vel_gauss[j] += v[j] - w[j];
	
      }
    diffuse_conductivity *= lumping_factor;//扩散系数平均
    //density *= lumping_factor;//密度
    //specific_heat *= lumping_factor;//比热
    concen_flux *= lumping_factor;//物质密度流
    proj *= lumping_factor;
    ms_vel_gauss *= lumping_factor;

    double c1 = 4.00;
    double c2 = 2.00;
    double h = sqrt(2.00 * Area);//单元尺寸？
    double norm_u = norm_2(ms_vel_gauss);//第二范数
    double tau1;
    if(Stationary==1){
      //tau1 = ( c1 * diffuse_conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);//密度*比热=1
      tau1 = ( c1 * diffuse_conductivity + c2 * (norm_u + 1e-6) * h);//可能为了稳定化求解
    }
    else{
      //tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * diffuse_conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);
      tau1 = (h * h) / (BDFcoeffs[0] * h * h + c1 * diffuse_conductivity + c2 * (norm_u + 1e-6) * h);
    }        
    
    double p1 = msDN_DX(0, 0) * GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(1, 0) * GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(2, 0) * GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar);
    double p2 = msDN_DX(0, 1) * GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(1, 1) * GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(2, 1) * GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar);
    grad_g[0] = p1;
    grad_g[1] = p2;//未知量的梯度
    
    
    double res = (inner_prod(ms_vel_gauss, grad_g)); //+ 0.333333333333333 * (t0media+t1media+t2media)*(1/dt)*density*diffuse_conductivity;
    double aux_res = fabs(res - proj);
    if (fabs(res) > aux_res)
      res = aux_res;
    //        res -= proj;
    //res *= density*specific_heat;//密度*比热=1
    double norm_grad = norm_2(grad_g);
    double k_aux = fabs(res) / (norm_grad + 1e-6);//残差与梯度范数的比值
    k_aux *= 0.707;
    
    noalias(First) = outer_prod(ms_vel_gauss, trans(ms_vel_gauss));
    First /= ((norm_u + 1e-6)*(norm_u + 1e-6));
    noalias(Second) = Identity - First;
    noalias(Third) = prod(Second, trans(msDN_DX));
    
    
    //CONCENTRATION CONTRIBUTION TO THE STIFFNESS MATRIX
    noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
    //noalias(rLeftHandSideMatrix) = (density * specific_heat) * outer_prod(msN, ms_u_DN);//密度*比热=1
    noalias(rLeftHandSideMatrix) = outer_prod(msN, ms_u_DN);//组装刚度阵
    
    //CONCENTRATION STABILIZING CONTRIBUTION (Suu)
    //noalias(rLeftHandSideMatrix) += density * specific_heat * density * specific_heat * tau1 * outer_prod(ms_u_DN, ms_u_DN);//密度*比热=1
    noalias(rLeftHandSideMatrix) += tau1 * outer_prod(ms_u_DN, ms_u_DN);//组装刚度阵
    
    //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX//粘性
    noalias(rLeftHandSideMatrix) += (diffuse_conductivity * prod(msDN_DX, trans(msDN_DX)) + k_aux * h * prod(msDN_DX, Third));
    
    //INERTIA CONTRIBUTION
    if(Stationary!=1){//非稳态问题
      //noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * (density * specific_heat) * msMassFactors;
      noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMassFactors;
    }
    // RHS = Fext (concentration per unit mass)//内热源
    //noalias(rRightHandSideVector) = (concen_flux * density) * msN;
    noalias(rRightHandSideVector) = (concen_flux) * msN;//右手边向量

    
    //RHS += Suy * proj[component]
    //noalias(rRightHandSideVector) += density * specific_heat * density * specific_heat * (tau1 * proj) * ms_u_DN;//密度*比热=1
    noalias(rRightHandSideVector) += (tau1 * proj) * ms_u_DN;//右手边向量
    
    //adding the inertia terms
    // RHS += M*vhistory
    //calculating the historical velocity
    
    if(Stationary!=1){//非稳态问题
      for (unsigned int iii = 0; iii < number_of_points; iii++)
	ms_concen_vec_np[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);//BDFcoeffs向后差分时间步长系数
      for (unsigned int step = 2; step < BDFcoeffs.size(); step++)
        {
	  for (unsigned int iii = 0; iii < number_of_points; iii++)
	    ms_concen_vec_np[iii] += BDFcoeffs[step] * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, step);
        }
      //noalias(rRightHandSideVector) -= prod(msMassFactors, ms_concen_vec_np * density * specific_heat);
      noalias(rRightHandSideVector) -= prod(msMassFactors, ms_concen_vec_np);//密度*比热=1
    }
    //subtracting the dirichlet term
    // RHS -= LHS*concentration
    for (unsigned int iii = 0; iii < number_of_points; iii++)
      ms_concen_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_concen_vec_np);
    
    
    rRightHandSideVector *= Area;
    rLeftHandSideMatrix *= Area;
    
    KRATOS_CATCH("");
  }
  
  

  //************************************************************************************
  //************************************************************************************

    void ConcenDiff2D::CalculateRightHandSide(VectorType& rRightHandSideVector, 
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void ConcenDiff2D::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
                int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        BoundedMatrix<double, 3, 2 > msDN_DX;
        array_1d<double, 3 > msN;
        array_1d<double, 2 > ms_vel_gauss;
        array_1d<double, 3 > ms_concen_vec_np;
        array_1d<double, 3 > ms_u_DN;

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
        ConcentrationDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
        const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
	const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();

        if (FractionalStepNumber == 2) //calculation of temperature convective projection
        {
            const unsigned int number_of_points = GetGeometry().size();
            const double lumping_factor = 1.00 / double(number_of_points);
            unsigned int TDim = 2;

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
                //const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
                const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
                for (unsigned int j = 0; j < TDim; j++)
                    ms_vel_gauss[j] += v[j] - w[j];

            }
            ms_vel_gauss *= lumping_factor;

            //calculating convective auxiliary vector
            noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
            double temp_conv = inner_prod(ms_u_DN, ms_concen_vec_np);
            temp_conv *= Area;

            for (unsigned int i = 0; i < number_of_points; i++)
            {
                GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
		GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable) += lumping_factor*temp_conv;
            }
        }
        KRATOS_CATCH("");
    }


    //************************************************************************************
    //************************************************************************************

    void ConcenDiff2D::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const 
    {
        ConcentrationDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONCENTRATION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    }


    void ConcenDiff2D::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const
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


