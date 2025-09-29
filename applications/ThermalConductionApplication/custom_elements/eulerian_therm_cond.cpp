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
#include "custom_elements/eulerian_therm_cond.h"
#include "thermal_conduction_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
//#include "applications/PhaseTransformationApplication/phase_transformation_application_variables.h"

namespace Kratos
{

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_TRY

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
        }

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::GetDofList(DofsVectorType &ElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_TRY

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
        }

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::CalculateLocalSystem(Matrix &rLeftHandSideMatrix,
                                                                                 Vector &rRightHandSideVector,
                                                                                 const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Resize of the Left and Right Hand side
        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); // false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); // false says not to preserve existing storage!!

        // Element variables
        ElementVariables Variables;
        this->InitializeEulerianElement(Variables, rCurrentProcessInfo);

        // Compute the geometry
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        array_1d<double, TNumNodes> N;
        double Volume;
        this->CalculateGeometry(DN_DX, Volume);//取一个中心导数，乘以体积近似积分

        // Getting the values of shape functions on Integration Points
        BoundedMatrix<double, TNumNodes, TNumNodes> Ncontainer;
        const GeometryType &Geom = this->GetGeometry();
        Ncontainer = Geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        // Getting the values of Current Process Info and computing the value of h
        this->GetNodalValues(Variables, rCurrentProcessInfo);

        // Some auxilary definitions
        BoundedMatrix<double, TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); // terms multiplying dphi/dt
        bounded_matrix<double, TNumNodes, TDim> tmp;

        // Gauss points and Number of nodes coincides in this case.
        for(unsigned int igauss=0; igauss<TNumNodes; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);

            //terms multiplying dphi/dt (aux1)
            noalias(aux1) += outer_prod(N, N);//节点值近似 Gauss 点，用平均法求积分

        }
        //adding the second and third term in the formulation
        noalias(rLeftHandSideMatrix)  = (Variables.dt_inv*Variables.density*Variables.specific_heat)*aux1;
        noalias(rRightHandSideVector) = (Variables.dt_inv*Variables.density*Variables.specific_heat)*prod(aux1,Variables.t_old);
        
        // adding the diffusion
        noalias(rLeftHandSideMatrix) += (Variables.conductivity * prod(DN_DX, trans(DN_DX))) * static_cast<double>(TNumNodes);

        // phase transformation latent, treat like volumn source (affecting the RHS only)
        const array_1d<double, TNumNodes> phase_latent = Variables.austenize_latent + Variables.eutectoid_latent + Variables.bainite_latent + Variables.martensite_latent;

        // volume source terms (affecting the RHS only)
        noalias(rRightHandSideVector) += prod(aux1, Variables.volumetric_source + phase_latent);

        // take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Variables.t_now);

        rRightHandSideVector *= Volume / static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume / static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Eulerian ThermCond Element")
    }

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::InitializeEulerianElement(ElementVariables &rVariables, const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA]; // Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        rVariables.dt_inv = 1.0 / delta_t;
        rVariables.lumping_factor = 1.00 / double(TNumNodes);

        rVariables.conductivity = 0.0;
        rVariables.specific_heat = 0.0;
        rVariables.density = 0.0;

        // 多相材料的质量分数初始化，设为0.0
        // rVariables.mass_fraction_m1 = 0.0;
        // rVariables.mass_fraction_m2 = 0.0;
        // rVariables.mass_fraction_m3 = 0.0;
        // rVariables.mass_fraction_m4 = 0.0;
        // rVariables.mass_fraction_m5 = 0.0;
        // rVariables.mass_fraction_m6 = 0.0;
        // rVariables.mass_fraction_m7 = 0.0;

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::CalculateGeometry(BoundedMatrix<double, TNumNodes, TDim> &rDN_DX, double &rVolume)
    {

        const GeometryType &Geom = this->GetGeometry();

        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType &integration_points = Geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1);
        const unsigned int NumGPoints = integration_points.size();
        rVolume = Geom.Area();
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, GeometryData::IntegrationMethod::GI_GAUSS_1);

        noalias(rDN_DX) = DN_DXContainer[0];
    }

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::GetNodalValues(ElementVariables &rVariables, const ProcessInfo &rCurrentProcessInfo) const
    {
        const Properties &r_material_properties = this->GetProperties();
        const auto& r_geom = this->GetGeometry();
        
        // check multiphase
        bool IsDefinedMultiPhase = false;
        if(r_material_properties.Has(IS_DEFINED_MULTIPHASE))
            IsDefinedMultiPhase = r_material_properties[IS_DEFINED_MULTIPHASE];

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rVariables.t_now[i] = r_geom[i].FastGetSolutionStepValue(TEMPERATURE);
            rVariables.t_old[i] = r_geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);

            // const auto &r_geom = GetGeometry();
            // Vector N(TNumNodes, 0.0);
            // N[i] = 1.0;
            // const auto r_N = N;
            // const auto &r_N = r_geom.ShapeFunctionsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 );//积分方法阶数

            rVariables.mass_fraction_current[i] = ZeroVector(7);
            rVariables.mass_fraction_previous[i] = ZeroVector(7);

            rVariables.volumetric_source[i] = 0.0;
            rVariables.austenize_latent[i] = 0.0;
            rVariables.eutectoid_latent[i] = 0.0;
            rVariables.bainite_latent[i] = 0.0;
            rVariables.martensite_latent[i] = 0.0;

            // 多相分数
            if (IsDefinedMultiPhase)
            {
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

            //密度
            if (IsDefinedMultiPhase)
            {
                double density_m1 = r_material_properties.GetValue(DENSITY_M1, r_geom[i], rCurrentProcessInfo);
                double density_m2 = r_material_properties.GetValue(DENSITY_M2, r_geom[i], rCurrentProcessInfo);
                double density_m3 = r_material_properties.GetValue(DENSITY_M3, r_geom[i], rCurrentProcessInfo);
                double density_m4 = r_material_properties.GetValue(DENSITY_M4, r_geom[i], rCurrentProcessInfo);
                double density_m5 = r_material_properties.GetValue(DENSITY_M5, r_geom[i], rCurrentProcessInfo);
                double density_m6 = r_material_properties.GetValue(DENSITY_M6, r_geom[i], rCurrentProcessInfo);
                double density_m7 = r_material_properties.GetValue(DENSITY_M7, r_geom[i], rCurrentProcessInfo);
                // std::cout << GetGeometry()[i].FastGetSolutionStepValue(my_settings->GetDensityM1Variable())<< std::endl;

                rVariables.density += (density_m1 * rVariables.mass_fraction_current[i][1] +
                                       density_m2 * rVariables.mass_fraction_current[i][2] +
                                       density_m3 * rVariables.mass_fraction_current[i][3] +
                                       density_m4 * rVariables.mass_fraction_current[i][4] +
                                       density_m5 * rVariables.mass_fraction_current[i][5] +
                                       density_m6 * rVariables.mass_fraction_current[i][6] +
                                       density_m7 * rVariables.mass_fraction_current[i][7]);
                // rVariables.density += (density_m1_sum + density_m2_sum + density_m3_sum + density_m4_sum + density_m5_sum + density_m6_sum + density_m7_sum);
            }else{
                if(r_material_properties.HasAccessor(DENSITY)){
                    double density_node = r_material_properties.GetValue(DENSITY, r_geom[i], rCurrentProcessInfo);
                    rVariables.density += density_node; // 密度计算
                }else{
                    rVariables.density += r_geom[i].FastGetSolutionStepValue(DENSITY);
                }
            }
            if (rVariables.density == 0.0){
                rVariables.density += 1.0;
                KRATOS_WARNING(this->Info()) << this->Id() << ":: No Density variable obtained or is 0.0, take default as 1.0" << std::endl;
            }

            // 比热
            if (IsDefinedMultiPhase)
            {
                double specific_heat_m1 = r_material_properties.GetValue(SPECIFIC_HEAT_M1, r_geom[i], rCurrentProcessInfo);
                double specific_heat_m2 = r_material_properties.GetValue(SPECIFIC_HEAT_M2, r_geom[i], rCurrentProcessInfo);
                double specific_heat_m3 = r_material_properties.GetValue(SPECIFIC_HEAT_M3, r_geom[i], rCurrentProcessInfo);
                double specific_heat_m4 = r_material_properties.GetValue(SPECIFIC_HEAT_M4, r_geom[i], rCurrentProcessInfo);
                double specific_heat_m5 = r_material_properties.GetValue(SPECIFIC_HEAT_M5, r_geom[i], rCurrentProcessInfo);
                double specific_heat_m6 = r_material_properties.GetValue(SPECIFIC_HEAT_M6, r_geom[i], rCurrentProcessInfo);
                double specific_heat_m7 = r_material_properties.GetValue(SPECIFIC_HEAT_M7, r_geom[i], rCurrentProcessInfo);

                rVariables.specific_heat += (specific_heat_m1 * rVariables.mass_fraction_current[i][1]+
                                                specific_heat_m2 * rVariables.mass_fraction_current[i][2]+
                                                specific_heat_m3 * rVariables.mass_fraction_current[i][3]+
                                                specific_heat_m4 * rVariables.mass_fraction_current[i][4]+
                                                specific_heat_m5 * rVariables.mass_fraction_current[i][5]+
                                                specific_heat_m6 * rVariables.mass_fraction_current[i][6]+
                                                specific_heat_m7 * rVariables.mass_fraction_current[i][7]);
                // rVariables.specific_heat += (specific_heat_m1_sum + specific_heat_m2_sum + specific_heat_m3_sum + specific_heat_m4_sum + specific_heat_m5_sum + specific_heat_m6_sum + specific_heat_m7_sum);
            }else{
                if(r_material_properties.HasAccessor(SPECIFIC_HEAT)){
                    double specific_heat_node = r_material_properties.GetValue(SPECIFIC_HEAT, r_geom[i], rCurrentProcessInfo);
                    rVariables.specific_heat += specific_heat_node; // 比热计算
                }else{
                    rVariables.specific_heat += r_geom[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
                }
            }
            if(rVariables.specific_heat == 0.0){
                rVariables.specific_heat += 1.0;
                KRATOS_WARNING(this->Info()) << this->Id() << ":: No Specific Heat variable obtained or is 0.0, take default as 1.0" << std::endl;
            }

            if (IsDefinedMultiPhase)
            {
                double conductivity_m1 = r_material_properties.GetValue(CONDUCTIVITY_M1, r_geom[i], rCurrentProcessInfo);
                double conductivity_m2 = r_material_properties.GetValue(CONDUCTIVITY_M2, r_geom[i], rCurrentProcessInfo);
                double conductivity_m3 = r_material_properties.GetValue(CONDUCTIVITY_M3, r_geom[i], rCurrentProcessInfo);
                double conductivity_m4 = r_material_properties.GetValue(CONDUCTIVITY_M4, r_geom[i], rCurrentProcessInfo);
                double conductivity_m5 = r_material_properties.GetValue(CONDUCTIVITY_M5, r_geom[i], rCurrentProcessInfo);
                double conductivity_m6 = r_material_properties.GetValue(CONDUCTIVITY_M6, r_geom[i], rCurrentProcessInfo);
                double conductivity_m7 = r_material_properties.GetValue(CONDUCTIVITY_M7, r_geom[i], rCurrentProcessInfo);

                rVariables.conductivity += (conductivity_m1 * rVariables.mass_fraction_current[i][1] +
                                            conductivity_m2 * rVariables.mass_fraction_current[i][2] +
                                            conductivity_m3 * rVariables.mass_fraction_current[i][3] +
                                            conductivity_m4 * rVariables.mass_fraction_current[i][4] +
                                            conductivity_m5 * rVariables.mass_fraction_current[i][5] +
                                            conductivity_m6 * rVariables.mass_fraction_current[i][6] +
                                            conductivity_m7 * rVariables.mass_fraction_current[i][7]);

                // rVariables.conductivity += (conductivity_m1_sum + conductivity_m2_sum + conductivity_m3_sum + conductivity_m4_sum + conductivity_m5_sum + conductivity_m6_sum + conductivity_m7_sum);
            }
            else
            {
                if(r_material_properties.HasAccessor(CONDUCTIVITY)){
                    double conductivity_node = r_material_properties.GetValue(CONDUCTIVITY, r_geom[i], rCurrentProcessInfo);
                    rVariables.conductivity += conductivity_node; // 导热系数计算
                }else{
                    rVariables.conductivity += r_geom[i].FastGetSolutionStepValue(CONDUCTIVITY);
                }
            }
            // if not, then the conductivity = 0
            if(rVariables.conductivity == 0.0){
                KRATOS_WARNING(this->Info()) << this->Id() << ":: No Conductivity variable obtained or is 0.0, take default as 1.0" << std::endl;
            }

            //体积热流
            rVariables.volumetric_source[i] += r_geom[i].FastGetSolutionStepValue(HEAT_FLUX);
            
            // 得到相变潜热值，###还需要根据相变类型确定相变分数增量!! ###
            // if(rCurrentProcessInfo[IS_PHASE_LATENT_CALCULATED] == true){
            if(rCurrentProcessInfo[IS_DEFINED_MULTIPHASE] == true){
            
                double austenize_latent_node = r_material_properties.GetValue(LATENT_AUSTENIZE, r_geom[i], rCurrentProcessInfo);
                rVariables.austenize_latent[i] += austenize_latent_node*(rVariables.mass_fraction_current[i][1] - rVariables.mass_fraction_previous[i][1]);

                double eutectoid_latent_node = r_material_properties.GetValue(LATENT_EUTECTOID_DECOMPOSITION, r_geom[i], rCurrentProcessInfo);
                rVariables.eutectoid_latent[i] += eutectoid_latent_node*(rVariables.mass_fraction_current[i][2] - rVariables.mass_fraction_previous[i][2]);

                double bainite_latent_node = r_material_properties.GetValue(LATENT_BAINITE_TRANSFORMATION, r_geom[i], rCurrentProcessInfo);
                rVariables.bainite_latent[i] += bainite_latent_node*(rVariables.mass_fraction_current[i][3] - rVariables.mass_fraction_previous[i][3]);

                double martensite_latent_node = r_material_properties.GetValue(LATENT_MARTENSITE_TRANSFORMATION, r_geom[i], rCurrentProcessInfo);
                rVariables.martensite_latent[i] += martensite_latent_node*(rVariables.mass_fraction_current[i][4] - rVariables.mass_fraction_previous[i][4]);
            }
            
        }

        rVariables.conductivity *= rVariables.lumping_factor;
        rVariables.density *= rVariables.lumping_factor;
        rVariables.specific_heat *= rVariables.lumping_factor;
    }

    int UpdateMassFraction(Kratos::PhaseTransformationLaw::Parameters &rValues, double *mfinc){
        const auto& rMaterialProperties = rValues.GetMaterialProperties();
        const auto& rElementGeometry = rValues.GetElementGeometry();
        const auto& rCurrentProcessInfo = rValues.GetProcessInfo();

        // Update the mass fraction of each phase
        // if(rCurrentProcessInfo[IS_PHASE_TRANSITION_CALCULATED] == true){
            
        // }

        return 0;
    }

    //----------------------------------------------------------------------------------------

    template <unsigned int TDim, unsigned int TNumNodes>
    void EulerianThermalConductionElement<TDim, TNumNodes>::CalculateRightHandSide(
        VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
    {
        Matrix LeftHandSide;
        this->CalculateLocalSystem(LeftHandSide, rRightHandSideVector, rCurrentProcessInfo);
    }

    //----------------------------------------------------------------------------------------

    template class EulerianThermalConductionElement<2, 3>;
    template class EulerianThermalConductionElement<2, 4>;
    template class EulerianThermalConductionElement<3, 4>;
    template class EulerianThermalConductionElement<3, 8>;

}
