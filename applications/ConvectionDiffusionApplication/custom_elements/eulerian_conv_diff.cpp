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
#include "custom_elements/eulerian_conv_diff.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::CalculateLocalSystem(Matrix& rLeftHandSideMatrix,
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
        double h = this->ComputeH(DN_DX);

        //Computing the divergence
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for(unsigned int k=0; k<TDim; k++)
            {
                Variables.div_v += DN_DX(i,k)*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
        }

        //Some auxilary definitions
        BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
        BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
        bounded_matrix<double,TNumNodes, TDim> tmp;

        // Gauss points and Number of nodes coincides in this case.
        for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //obtain the velocity in the middle of the tiem step
            array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                 for(unsigned int k=0; k<TDim; k++)
                    vel_gauss[k] += N[i]*(Variables.v[i][k]*Variables.theta + Variables.vold[i][k]*(1.0-Variables.theta));
            }
            const double norm_vel = norm_2(vel_gauss);
            array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

            const double tau = this->CalculateTau(Variables,norm_vel,h);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, N);
            noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

            //terms which multiply the gradient of phi
            noalias(aux2) += (1.0+tau*Variables.beta*Variables.div_v)*outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);
        }

        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (Variables.dt_inv*Variables.density*Variables.specific_heat + Variables.theta*Variables.beta*Variables.div_v)*aux1;
        noalias(rRightHandSideVector) = (Variables.dt_inv*Variables.density*Variables.specific_heat - (1.0-Variables.theta)*Variables.beta*Variables.div_v)*prod(aux1,Variables.phi_old);

        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (Variables.conductivity * Variables.theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((Variables.conductivity * (1.0-Variables.theta) * prod(DN_DX, trans(DN_DX))),Variables.phi_old)*static_cast<double>(TNumNodes) ;

        //terms in aux2
        noalias(rLeftHandSideMatrix) += Variables.density*Variables.specific_heat*Variables.theta*aux2;
        noalias(rRightHandSideVector) -= Variables.density*Variables.specific_heat*(1.0-Variables.theta)*prod(aux2,Variables.phi_old);

        // volume source terms (affecting the RHS only)
        noalias(rRightHandSideVector) += prod(aux1, Variables.volumetric_source);

        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Eulerian ConvDiff Element")
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA]; //Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        rVariables.dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
		rVariables.lumping_factor = 1.00 / double(TNumNodes);

        rVariables.conductivity = 0.0;
        rVariables.specific_heat = 0.0;
        rVariables.density = 0.0;
        rVariables.beta = 0.0;
        rVariables.div_v = 0.0;


        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::CalculateGeometry(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume)
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
    double EulerianConvectionDiffusionElement< TDim, TNumNodes >::ComputeH(BoundedMatrix<double,TNumNodes,TDim >& DN_DX)
    {
        double h=0.0;

        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement<TDim,TNumNodes>::GetNodalValues(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) const
    {
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

        //////storing locally the flags to avoid repeated check in the nodal loops
        const bool IsDefinedVelocityVariable = my_settings->IsDefinedVelocityVariable();
        const bool IsDefinedMeshVelocityVariable = my_settings->IsDefinedMeshVelocityVariable();
        const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
        const bool IsDefinedLatentAustenizeVariable = my_settings->IsDefinedLatentAustenizeVariable();
        const bool IsDefinedLatentEutectoidVariable = my_settings->IsDefinedLatentEutectoidVariable();
        const bool IsDefinedLatentBainiteVariable = my_settings->IsDefinedLatentBainiteVariable();
        const bool IsDefinedLatentMartensiteVariable = my_settings->IsDefinedLatentMartensiteVariable();
        const bool IsDefinedSpecificHeatVariableVariable = my_settings->IsDefinedSpecificHeatVariable();
        const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();
        const bool IsDefinedVolumeSourceVariable = my_settings->IsDefinedVolumeSourceVariable();
        //check multiphase
        //const bool IsDefinedMultiPhase = (my_settings->IsDefinedMassFractionM1Variable()) || (my_settings->IsDefinedMassFractionM2Variable()) || (my_settings->IsDefinedMassFractionM3Variable()) || (my_settings->IsDefinedMassFractionM4Variable()) || (my_settings->IsDefinedMassFractionM5Variable()) || (my_settings->IsDefinedMassFractionM6Variable()) || (my_settings->IsDefinedMassFractionM7Variable());
        const bool IsDefinedMultiPhase = false;
        
        const Properties& r_material_properties = this->GetProperties();
        const auto& r_geom = this->GetGeometry();

        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        const Variable<double>& rMassFractionM1Var = my_settings->GetMassFractionM1Variable();
        const Variable<double>& rMassFractionM2Var = my_settings->GetMassFractionM2Variable();
        const Variable<double>& rMassFractionM3Var = my_settings->GetMassFractionM3Variable();
        const Variable<double>& rMassFractionM4Var = my_settings->GetMassFractionM4Variable();
        const Variable<double>& rMassFractionM5Var = my_settings->GetMassFractionM5Variable();
        const Variable<double>& rMassFractionM6Var = my_settings->GetMassFractionM6Variable();
        const Variable<double>& rMassFractionM7Var = my_settings->GetMassFractionM7Variable();


        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.temperature_current[i] = r_geom[i].FastGetSolutionStepValue(rUnknownVar);
            rVariables.temperature_previous[i] = r_geom[i].FastGetSolutionStepValue(rUnknownVar,1);
            
            rVariables.phi[i] = r_geom[i].FastGetSolutionStepValue(rUnknownVar);
            rVariables.phi_old[i] = r_geom[i].FastGetSolutionStepValue(rUnknownVar,1);
            //dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

			rVariables.v[i]=ZeroVector(3);
			rVariables.vold[i]=ZeroVector(3);

            rVariables.mass_fraction_current[i]=ZeroVector(7);
	        rVariables.mass_fraction_previous[i]=ZeroVector(7);

            rVariables.volumetric_source[i] = 0.0;
            rVariables.austenize_latent[i] = 0.0;
            rVariables.eutectoid_latent[i] = 0.0;
            rVariables.bainite_latent[i] = 0.0;
            rVariables.martensite_latent[i] = 0.0;

            //多相分数
            if (IsDefinedMultiPhase){
                rVariables.mass_fraction_current[i][1] = r_geom[i].FastGetSolutionStepValue(rMassFractionM1Var);
                rVariables.mass_fraction_current[i][2] = r_geom[i].FastGetSolutionStepValue(rMassFractionM2Var);
                rVariables.mass_fraction_current[i][3] = r_geom[i].FastGetSolutionStepValue(rMassFractionM3Var);
                rVariables.mass_fraction_current[i][4] = r_geom[i].FastGetSolutionStepValue(rMassFractionM4Var);
                rVariables.mass_fraction_current[i][5] = r_geom[i].FastGetSolutionStepValue(rMassFractionM5Var);
                rVariables.mass_fraction_current[i][6] = r_geom[i].FastGetSolutionStepValue(rMassFractionM6Var);
                rVariables.mass_fraction_current[i][7] = r_geom[i].FastGetSolutionStepValue(rMassFractionM7Var);
                        
                rVariables.mass_fraction_previous[i][1] = r_geom[i].FastGetSolutionStepValue(rMassFractionM1Var, 1);
                rVariables.mass_fraction_previous[i][2] = r_geom[i].FastGetSolutionStepValue(rMassFractionM2Var, 1);
                rVariables.mass_fraction_previous[i][3] = r_geom[i].FastGetSolutionStepValue(rMassFractionM3Var, 1);
                rVariables.mass_fraction_previous[i][4] = r_geom[i].FastGetSolutionStepValue(rMassFractionM4Var, 1);
                rVariables.mass_fraction_previous[i][5] = r_geom[i].FastGetSolutionStepValue(rMassFractionM5Var, 1);
                rVariables.mass_fraction_previous[i][6] = r_geom[i].FastGetSolutionStepValue(rMassFractionM6Var, 1);
                rVariables.mass_fraction_previous[i][7] = r_geom[i].FastGetSolutionStepValue(rMassFractionM7Var, 1);
                  
            }

            if (IsDefinedVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
				  rVariables.v[i] = r_geom[i].FastGetSolutionStepValue(rVelocityVar);
				  rVariables.vold[i] = r_geom[i].FastGetSolutionStepValue(rVelocityVar,1);
				  //active_convection=true;
			}

			if (IsDefinedMeshVelocityVariable)
            {
				  const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
				  rVariables.v[i] -= r_geom[i].FastGetSolutionStepValue(rMeshVelocityVar);
				  rVariables.vold[i] -= r_geom[i].FastGetSolutionStepValue(rMeshVelocityVar,1);
				  //active_convection=true;
			}

			if (IsDefinedDensityVariable)
			{
				if (IsDefinedMultiPhase)
		        {
                    double density_m1 = r_material_properties.GetValue(my_settings->GetDensityM1Variable(), r_geom[i], rCurrentProcessInfo);
                    double density_m2 = r_material_properties.GetValue(my_settings->GetDensityM2Variable(), r_geom[i], rCurrentProcessInfo);
                    double density_m3 = r_material_properties.GetValue(my_settings->GetDensityM3Variable(), r_geom[i], rCurrentProcessInfo);
                    double density_m4 = r_material_properties.GetValue(my_settings->GetDensityM4Variable(), r_geom[i], rCurrentProcessInfo);
                    double density_m5 = r_material_properties.GetValue(my_settings->GetDensityM5Variable(), r_geom[i], rCurrentProcessInfo);
                    double density_m6 = r_material_properties.GetValue(my_settings->GetDensityM6Variable(), r_geom[i], rCurrentProcessInfo);
                    double density_m7 = r_material_properties.GetValue(my_settings->GetDensityM7Variable(), r_geom[i], rCurrentProcessInfo);
                    //std::cout << r_geom[i].FastGetSolutionStepValue(my_settings->GetDensityM1Variable())<< std::endl;
                
                    rVariables.density += (density_m1*r_geom[i].FastGetSolutionStepValue(rMassFractionM1Var) + 
                                        density_m2*r_geom[i].FastGetSolutionStepValue(rMassFractionM2Var) + 
                                        density_m3*r_geom[i].FastGetSolutionStepValue(rMassFractionM3Var) + 
                                        density_m4*r_geom[i].FastGetSolutionStepValue(rMassFractionM4Var) + 
                                        density_m5*r_geom[i].FastGetSolutionStepValue(rMassFractionM5Var) + 
                                        density_m6*r_geom[i].FastGetSolutionStepValue(rMassFractionM6Var) + 
                                        density_m7*r_geom[i].FastGetSolutionStepValue(rMassFractionM7Var) );
                    //rVariables.density += (density_m1_sum + density_m2_sum + density_m3_sum + density_m4_sum + density_m5_sum + density_m6_sum + density_m7_sum);
			    }else{
                    // std::cout << "#### density" << std::endl;
                    // std::cout << "#### mat_properties" << r_material_properties << std::endl;
                    // std::cout << "#### varialbles" << my_settings->GetDensityVariable() << std::endl;
                    // std::cout << "#### geometry" << r_geom << std::endl;
                    // std::cout << "#### geometry[i]" << r_geom[i] << std::endl;
                    double density_node = r_material_properties.GetValue(my_settings->GetDensityVariable(), r_geom[i], rCurrentProcessInfo);
                    rVariables.density += density_node;//密度计算
                }
                if (rVariables.density == 0.0){
                    const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
				    rVariables.density += r_geom[i].FastGetSolutionStepValue(rDensityVar);
                }
			}
			else
				rVariables.density += 1.0;
            
            const bool first_computation = (rCurrentProcessInfo[NL_ITERATION_NUMBER] == 1 && rCurrentProcessInfo[STEP] == 1) ? true : false;
            if(first_computation){
                std::ofstream output_file("test_nonhistorical_historical.txt", std::ios::trunc);
                output_file.close();
            }
            
			if (IsDefinedSpecificHeatVariableVariable)
			{
				if (IsDefinedMultiPhase)
		        {
                    double specific_heat_m1 = r_material_properties.GetValue(my_settings->GetSpecificHeatM1Variable(), r_geom[i], rCurrentProcessInfo);
                    double specific_heat_m2 = r_material_properties.GetValue(my_settings->GetSpecificHeatM2Variable(), r_geom[i], rCurrentProcessInfo);
                    double specific_heat_m3 = r_material_properties.GetValue(my_settings->GetSpecificHeatM3Variable(), r_geom[i], rCurrentProcessInfo);
                    double specific_heat_m4 = r_material_properties.GetValue(my_settings->GetSpecificHeatM4Variable(), r_geom[i], rCurrentProcessInfo);
                    double specific_heat_m5 = r_material_properties.GetValue(my_settings->GetSpecificHeatM5Variable(), r_geom[i], rCurrentProcessInfo);
                    double specific_heat_m6 = r_material_properties.GetValue(my_settings->GetSpecificHeatM6Variable(), r_geom[i], rCurrentProcessInfo);
                    double specific_heat_m7 = r_material_properties.GetValue(my_settings->GetSpecificHeatM7Variable(), r_geom[i], rCurrentProcessInfo);
                    //std::cout << r_geom[i].FastGetSolutionStepValue(my_settings->GetDensityM1Variable())<< std::endl;
                
                    rVariables.specific_heat += (specific_heat_m1*r_geom[i].FastGetSolutionStepValue(rMassFractionM1Var) + 
				                                 specific_heat_m2*r_geom[i].FastGetSolutionStepValue(rMassFractionM2Var) + 
				                                 specific_heat_m3*r_geom[i].FastGetSolutionStepValue(rMassFractionM3Var) + 
				                                 specific_heat_m4*r_geom[i].FastGetSolutionStepValue(rMassFractionM4Var) + 
				                                 specific_heat_m5*r_geom[i].FastGetSolutionStepValue(rMassFractionM5Var) + 
				                                 specific_heat_m6*r_geom[i].FastGetSolutionStepValue(rMassFractionM6Var) + 
				                                 specific_heat_m7*r_geom[i].FastGetSolutionStepValue(rMassFractionM7Var) );
                    //rVariables.density += (density_m1_sum + density_m2_sum + density_m3_sum + density_m4_sum + density_m5_sum + density_m6_sum + density_m7_sum);
			    }else{
                    double specific_heat_node = r_material_properties.GetValue(my_settings->GetSpecificHeatVariable(), r_geom[i], rCurrentProcessInfo);
                    rVariables.specific_heat += specific_heat_node;//比热计算
                    if (r_geom[i].Id() == 1){
                        std::ofstream output_file;
                        output_file.open("test_nonhistorical_historical.txt", std::ios::app);
                        output_file<<"### STEP"<<rCurrentProcessInfo[STEP]<<"### ITERATION" <<rCurrentProcessInfo[NL_ITERATION_NUMBER] << std::endl;
                        output_file<<"###UnkownVar fast: "<<r_geom[i].FastGetSolutionStepValue(rUnknownVar)<<std::endl;
                        output_file<<"###UnkownVar Get: "<<r_geom[i].GetValue(rUnknownVar)<<std::endl;
                        output_file<<"###Temperature FastGet: "<<r_geom[i].FastGetSolutionStepValue(TEMPERATURE)<<std::endl;
                        output_file<<"###Temperature Get: "<<r_geom[i].GetValue(TEMPERATURE)<<std::endl;
                        output_file<<"###SpecificHeat accessor type: "<<r_material_properties.GetAccessor(SPECIFIC_HEAT).Info()<<std::endl;
                        output_file<<"###SpecificHeat: "<<specific_heat_node<<std::endl;
                        output_file<<"###SpecificHeat fast: "<<r_geom[i].FastGetSolutionStepValue(my_settings->GetSpecificHeatVariable())<<std::endl;
                        output_file<<"###SpecificHeat fast SPECIFIC_HEAT: "<<r_geom[i].FastGetSolutionStepValue(SPECIFIC_HEAT)<<std::endl;
                        output_file<<"###SpecificHeat Get: "<<r_geom[i].GetValue(my_settings->GetSpecificHeatVariable())<<std::endl;
                        output_file<<"###SpecificHeat Get SPECIFIC_HEAT: "<<r_geom[i].GetValue(SPECIFIC_HEAT)<<std::endl;
                        output_file.close();
                    }
                }
                if (rVariables.specific_heat == 0.0){
                    const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
				    rVariables.specific_heat += r_geom[i].FastGetSolutionStepValue(rSpecificHeatVar);
                }
			}
			else
				rVariables.specific_heat += 1.0;

			if (IsDefinedDiffusionVariable)
			{
				if (IsDefinedMultiPhase)
		        {
                    double conductivity_m1 = r_material_properties.GetValue(my_settings->GetDiffusionM1Variable(), r_geom[i], rCurrentProcessInfo);
                    double conductivity_m2 = r_material_properties.GetValue(my_settings->GetDiffusionM2Variable(), r_geom[i], rCurrentProcessInfo);
                    double conductivity_m3 = r_material_properties.GetValue(my_settings->GetDiffusionM3Variable(), r_geom[i], rCurrentProcessInfo);
                    double conductivity_m4 = r_material_properties.GetValue(my_settings->GetDiffusionM4Variable(), r_geom[i], rCurrentProcessInfo);
                    double conductivity_m5 = r_material_properties.GetValue(my_settings->GetDiffusionM5Variable(), r_geom[i], rCurrentProcessInfo);
                    double conductivity_m6 = r_material_properties.GetValue(my_settings->GetDiffusionM6Variable(), r_geom[i], rCurrentProcessInfo);
                    double conductivity_m7 = r_material_properties.GetValue(my_settings->GetDiffusionM7Variable(), r_geom[i], rCurrentProcessInfo);
                    //std::cout << r_geom[i].FastGetSolutionStepValue(my_settings->GetDensityM1Variable())<< std::endl;
                
                    rVariables.conductivity += (conductivity_m1*r_geom[i].FastGetSolutionStepValue(rMassFractionM1Var) + 
				                                conductivity_m2*r_geom[i].FastGetSolutionStepValue(rMassFractionM2Var) + 
				                                conductivity_m3*r_geom[i].FastGetSolutionStepValue(rMassFractionM3Var) + 
				                                conductivity_m4*r_geom[i].FastGetSolutionStepValue(rMassFractionM4Var) + 
				                                conductivity_m5*r_geom[i].FastGetSolutionStepValue(rMassFractionM5Var) + 
				                                conductivity_m6*r_geom[i].FastGetSolutionStepValue(rMassFractionM6Var) + 
				                                conductivity_m7*r_geom[i].FastGetSolutionStepValue(rMassFractionM7Var) );
                    //rVariables.density += (density_m1_sum + density_m2_sum + density_m3_sum + density_m4_sum + density_m5_sum + density_m6_sum + density_m7_sum);
			    }else{
                    double conductivity_node = r_material_properties.GetValue(my_settings->GetDiffusionVariable(), r_geom[i], rCurrentProcessInfo);
                    rVariables.conductivity += conductivity_node;//热导计算
                    if (r_geom[i].Id() == 1){
                        std::ofstream output_file;
                        output_file.open("test_nonhistorical_historical.txt", std::ios::app);
                        output_file<<"### STEP"<<rCurrentProcessInfo[STEP]<<"### ITERATION" <<rCurrentProcessInfo[NL_ITERATION_NUMBER] << std::endl;
                        output_file<<"###UnkownVar fast: "<<r_geom[i].FastGetSolutionStepValue(rUnknownVar)<<std::endl;
                        output_file<<"###UnkownVar Get: "<<r_geom[i].GetValue(rUnknownVar)<<std::endl;
                        output_file<<"###Temperature FastGet: "<<r_geom[i].FastGetSolutionStepValue(TEMPERATURE)<<std::endl;
                        output_file<<"###Temperature Get: "<<r_geom[i].GetValue(TEMPERATURE)<<std::endl;
                        output_file<<"###Conductivity: "<<conductivity_node<<std::endl;
                        output_file<<"###Conductivity fast: "<<r_geom[i].FastGetSolutionStepValue(my_settings->GetDiffusionVariable())<<std::endl;
                        output_file<<"###Conductivity fast SPECIFIC_HEAT: "<<r_geom[i].FastGetSolutionStepValue(CONDUCTIVITY)<<std::endl;
                        output_file<<"###Conductivity Get: "<<r_geom[i].GetValue(my_settings->GetDiffusionVariable())<<std::endl;
                        output_file<<"###Conductivity Get SPECIFIC_HEAT: "<<r_geom[i].GetValue(CONDUCTIVITY)<<std::endl;
                        output_file.close();
                    }
                }
                if (rVariables.conductivity == 0.0){
                    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
				    rVariables.conductivity += r_geom[i].FastGetSolutionStepValue(rDiffusionVar);
                }
			}
			//if not, then the conductivity = 0

            if (IsDefinedVolumeSourceVariable)
            {
                const Variable<double>& rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
                rVariables.volumetric_source[i] += r_geom[i].FastGetSolutionStepValue(rVolumeSourceVar);
            }
            //得到相变潜热值，###还需要根据相变类型确定相变分数增量!! ###
            if (IsDefinedLatentAustenizeVariable)
            {
                const Variable<double>& rLatentAustenizeVar = my_settings->GetLatentAustenizeVariable();
                rVariables.austenize_latent[i] += r_material_properties.GetValue(rLatentAustenizeVar, r_geom[i], rCurrentProcessInfo);
                rVariables.austenize_latent[i] *= ( rVariables.mass_fraction_current[i][1] - rVariables.mass_fraction_previous[i][1] );
            }
            if (IsDefinedLatentEutectoidVariable)
            {
                const Variable<double>& rLatentEutectoidVar = my_settings->GetLatentEutectoidVariable();
                rVariables.eutectoid_latent[i] += r_material_properties.GetValue(rLatentEutectoidVar, r_geom[i], rCurrentProcessInfo);
                rVariables.eutectoid_latent[i] *= ( rVariables.mass_fraction_current[i][2] - rVariables.mass_fraction_previous[i][2] );
            }
            if (IsDefinedLatentBainiteVariable)
            {
                const Variable<double>& rLatentBainiteVar = my_settings->GetLatentBainiteVariable();
                rVariables.bainite_latent[i] += r_material_properties.GetValue(rLatentBainiteVar, r_geom[i], rCurrentProcessInfo);
                rVariables.bainite_latent[i] *= ( rVariables.mass_fraction_current[i][3] - rVariables.mass_fraction_previous[i][2] );
            }
            if (IsDefinedLatentMartensiteVariable)
            {
                const Variable<double>& rLatentMartensiteVar = my_settings->GetLatentMartensiteVariable();
                rVariables.martensite_latent[i] += r_material_properties.GetValue(rLatentMartensiteVar, r_geom[i], rCurrentProcessInfo);
                rVariables.martensite_latent[i] *= ( rVariables.mass_fraction_current[i][4] - rVariables.mass_fraction_previous[i][4] );
            }
        }

        //array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
        //const double norm_grad = norm_2(grad_phi_halfstep);

        rVariables.conductivity *= rVariables.lumping_factor;
        rVariables.density *= rVariables.lumping_factor;
        rVariables.specific_heat *= rVariables.lumping_factor;
    }

    template< unsigned int TDim, unsigned int TNumNodes >
    double EulerianConvectionDiffusionElement<TDim,TNumNodes>::CalculateTau(const ElementVariables& rVariables, double norm_vel, double h)
    {
        // Dynamic part
        double inv_tau = rVariables.dyn_st_beta * rVariables.dt_inv;

        // Convection
        inv_tau += 2.0 * norm_vel / h + rVariables.beta*rVariables.div_v;

        // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
        inv_tau *= rVariables.density * rVariables.specific_heat;

        // Diffusion
        inv_tau += 4.0 * rVariables.conductivity / (h*h);

        // Limiting
        inv_tau = std::max(inv_tau, 1e-2);

        return (rVariables.density*rVariables.specific_heat) / inv_tau;
    }

//----------------------------------------------------------------------------------------

    template< unsigned int TDim, unsigned int TNumNodes >
    void EulerianConvectionDiffusionElement< TDim, TNumNodes >::CalculateRightHandSide(
        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        Matrix LeftHandSide;
        this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
    }

//----------------------------------------------------------------------------------------


template class EulerianConvectionDiffusionElement<2,3>;
template class EulerianConvectionDiffusionElement<2,4>;
template class EulerianConvectionDiffusionElement<3,4>;
template class EulerianConvectionDiffusionElement<3,8>;

}
