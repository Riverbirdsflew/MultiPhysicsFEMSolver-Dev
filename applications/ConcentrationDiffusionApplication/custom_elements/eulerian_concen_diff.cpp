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
#include "custom_elements/eulerian_concen_diff.h"
#include "concentration_diffusion_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement< TDim, TNumNodes >::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(CONCENTRATION).EquationId();
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement< TDim, TNumNodes >::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(CONCENTRATION);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement<TDim,TNumNodes>::CalculateLocalSystem(Matrix& rLeftHandSideMatrix,
                        Vector& rRightHandSideVector,
                        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize of the Left and Right Hand side
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!

        //Element variables
        ElementVariables Variables;
        this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

        // Compute the geometry
        BoundedMatrix<double,TNumNodes, TDim> DN_DX;
        array_1d<double,TNumNodes > N;
        double Volume;
        this-> CalculateGeometry(DN_DX,Volume);

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        const GeometryType& Geom = this->GetGeometry();
        Ncontainer = Geom.ShapeFunctionsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 );

        // Getting the values of Current Process Info and computing the value of h
        this-> GetNodalValues(Variables,rCurrentProcessInfo);

        //Some auxilary definitions
        BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dc/dt
        bounded_matrix<double,TNumNodes, TDim> tmp;

        // Gauss points and Number of nodes coincides in this case.
        for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //terms multiplying dc/dt (aux1)
            noalias(aux1) += outer_prod(N, N);
        }
        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (Variables.dt_inv)*aux1;
        noalias(rRightHandSideVector) = (Variables.dt_inv)*prod(aux1,Variables.c_old);

        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (Variables.diffuse_conductivity * Variables.theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);

        // volume source terms (affecting the RHS only)
        noalias(rRightHandSideVector) += prod(aux1, Variables.volumetric_source);

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.c_now);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Eulerian ConcenDiff Element")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement<TDim,TNumNodes>::InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA]; //Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
		rVariables.lumping_factor = 1.00 / double(TNumNodes);

        rVariables.diffuse_conductivity = 0.0;
        //密度 比热
        //rVariables.specific_heat = 0.0;
        //rVariables.density = 0.0;


        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement<TDim,TNumNodes>::CalculateGeometry(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume)
    {

        const GeometryType& Geom = this->GetGeometry();

        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( GeometryData::IntegrationMethod::GI_GAUSS_1 );
        const unsigned int NumGPoints = integration_points.size();
        rVolume = Geom.Area();
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,GeometryData::IntegrationMethod::GI_GAUSS_1);

        noalias( rDN_DX ) = DN_DXContainer[0];

    }


//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement<TDim,TNumNodes>::GetNodalValues(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) const
    {
        const Properties& r_material_properties = GetProperties();
        const auto &r_geom = GetGeometry();

        //check multiphase
        bool IsDefinedMultiPhase = false;
        if (r_material_properties.Has(IS_DEFINED_MULTIPHASE))
        {
            IsDefinedMultiPhase = r_material_properties[IS_DEFINED_MULTIPHASE];
        }

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            
            //const auto &r_N = r_geom.ShapeFunctionsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 );//积分方法阶数
            
            rVariables.c_now[i] = GetGeometry()[i].FastGetSolutionStepValue(CONCENTRATION);
            rVariables.c_old[i] = GetGeometry()[i].FastGetSolutionStepValue(CONCENTRATION,1);
            //dc_dt[i] = dt_inv*(c[i] - c_old [i];

            rVariables.volumetric_source[i] = 0.0;

            if (IsDefinedMultiPhase){
                rVariables.mass_fraction_current[i][1] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M1);
                rVariables.mass_fraction_current[i][2] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M2);
                rVariables.mass_fraction_current[i][3] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M3);
                rVariables.mass_fraction_current[i][4] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M4);
                rVariables.mass_fraction_current[i][5] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M5);
                rVariables.mass_fraction_current[i][6] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M6);
                rVariables.mass_fraction_current[i][7] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M7);

                rVariables.mass_fraction_previous[i][1] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M1, 1);
                rVariables.mass_fraction_previous[i][2] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M2, 1);
                rVariables.mass_fraction_previous[i][3] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M3, 1);
                rVariables.mass_fraction_previous[i][4] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M4, 1);
                rVariables.mass_fraction_previous[i][5] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M5, 1);
                rVariables.mass_fraction_previous[i][6] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M6, 1);
                rVariables.mass_fraction_previous[i][7] = r_geom[i].FastGetSolutionStepValue(MASS_FRACTION_M7, 1);
            }

            //扩散系数
            if (IsDefinedMultiPhase){
                double diffuse_conductivity_m1 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M1, r_geom[i], rCurrentProcessInfo);
                double diffuse_conductivity_m2 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M2, r_geom[i], rCurrentProcessInfo);
                double diffuse_conductivity_m3 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M3, r_geom[i], rCurrentProcessInfo);
                double diffuse_conductivity_m4 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M4, r_geom[i], rCurrentProcessInfo);
                double diffuse_conductivity_m5 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M5, r_geom[i], rCurrentProcessInfo);
                double diffuse_conductivity_m6 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M6, r_geom[i], rCurrentProcessInfo);
                double diffuse_conductivity_m7 = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY_M7, r_geom[i], rCurrentProcessInfo);
                
                rVariables.diffuse_conductivity += (diffuse_conductivity_m1 * rVariables.mass_fraction_current[i][1] +
                                                    diffuse_conductivity_m2 * rVariables.mass_fraction_current[i][2] +
                                                    diffuse_conductivity_m3 * rVariables.mass_fraction_current[i][3] +
                                                    diffuse_conductivity_m4 * rVariables.mass_fraction_current[i][4] +
                                                    diffuse_conductivity_m5 * rVariables.mass_fraction_current[i][5] +
                                                    diffuse_conductivity_m6 * rVariables.mass_fraction_current[i][6] +
                                                    diffuse_conductivity_m7 * rVariables.mass_fraction_current[i][7]);
                            
            }
            else{
                if(r_material_properties.HasAccessor(DIFFUSE_CONDUCTIVITY)){
                    double diffuse_conductivity_node = r_material_properties.GetValue(DIFFUSE_CONDUCTIVITY, r_geom[i], rCurrentProcessInfo);
                    rVariables.diffuse_conductivity += diffuse_conductivity_node; // 扩散系数计算
                }else{
                    rVariables.diffuse_conductivity += GetGeometry()[i].FastGetSolutionStepValue(DIFFUSE_CONDUCTIVITY);
                }
            }
            if(rVariables.diffuse_conductivity == 0.0)
            {
                KRATOS_WARNING(this->Info()) << this->Id() << ":: No Diffuse Conductivity variable obtained or is 0.0, take default as 1.0" << std::endl;
            }

			//if not, then the diffuse_conductivity = 0

            rVariables.volumetric_source[i] += GetGeometry()[i].FastGetSolutionStepValue(CONCENTRATION_FLUX);
        }

        rVariables.diffuse_conductivity *= rVariables.lumping_factor;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    int EulerianConcentrationDiffusionElement<TDim,TNumNodes>::UpdateCarbonContent(ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
    {
        return 0;
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConcentrationDiffusionElement< TDim, TNumNodes >::CalculateRightHandSide(
        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        Matrix LeftHandSide;
        this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
    }

//----------------------------------------------------------------------------------------


template class EulerianConcentrationDiffusionElement<2,3>;
template class EulerianConcentrationDiffusionElement<2,4>;
template class EulerianConcentrationDiffusionElement<3,4>;
template class EulerianConcentrationDiffusionElement<3,8>;

}
