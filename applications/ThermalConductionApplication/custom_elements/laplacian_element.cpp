// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
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
#include "custom_elements/laplacian_element.h"

#include "includes/checks.h"
#include "includes/define.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
LaplacianElement::LaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
LaplacianElement::LaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer LaplacianElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer LaplacianElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianElement>(NewId, pGeom, pProperties);
}

LaplacianElement::~LaplacianElement()
{
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_m1_var = r_settings.GetDiffusionM1Variable();
    const Variable<double>& r_diffusivity_m2_var = r_settings.GetDiffusionM2Variable();
    const Variable<double>& r_diffusivity_m3_var = r_settings.GetDiffusionM3Variable();
    const Variable<double>& r_diffusivity_m4_var = r_settings.GetDiffusionM4Variable();
    const Variable<double>& r_diffusivity_m5_var = r_settings.GetDiffusionM5Variable();
    const Variable<double>& r_diffusivity_m6_var = r_settings.GetDiffusionM6Variable();
    const Variable<double>& r_diffusivity_m7_var = r_settings.GetDiffusionM7Variable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();
    const Variable<double>& r_density_m1_var = r_settings.GetDensityM1Variable();
    const Variable<double>& r_density_m2_var = r_settings.GetDensityM2Variable();
    const Variable<double>& r_density_m3_var = r_settings.GetDensityM3Variable();
    const Variable<double>& r_density_m4_var = r_settings.GetDensityM4Variable();
    const Variable<double>& r_density_m5_var = r_settings.GetDensityM5Variable();
    const Variable<double>& r_density_m6_var = r_settings.GetDensityM6Variable();
    const Variable<double>& r_density_m7_var = r_settings.GetDensityM7Variable();
    const Variable<double>& r_specific_heat_m1_var = r_settings.GetSpecificHeatM1Variable();
    const Variable<double>& r_specific_heat_m2_var = r_settings.GetSpecificHeatM2Variable();
    const Variable<double>& r_specific_heat_m3_var = r_settings.GetSpecificHeatM3Variable();
    const Variable<double>& r_specific_heat_m4_var = r_settings.GetSpecificHeatM4Variable();
    const Variable<double>& r_specific_heat_m5_var = r_settings.GetSpecificHeatM5Variable();
    const Variable<double>& r_specific_heat_m6_var = r_settings.GetSpecificHeatM6Variable();
    const Variable<double>& r_specific_heat_m7_var = r_settings.GetSpecificHeatM7Variable();
    const Variable<double>& r_mass_fraction_m1_var = r_settings.GetMassFractionM1Variable();
    const Variable<double>& r_mass_fraction_m2_var = r_settings.GetMassFractionM2Variable();
    const Variable<double>& r_mass_fraction_m3_var = r_settings.GetMassFractionM3Variable();
    const Variable<double>& r_mass_fraction_m4_var = r_settings.GetMassFractionM4Variable();
    const Variable<double>& r_mass_fraction_m5_var = r_settings.GetMassFractionM5Variable();
    const Variable<double>& r_mass_fraction_m6_var = r_settings.GetMassFractionM6Variable();
    const Variable<double>& r_mass_fraction_m7_var = r_settings.GetMassFractionM7Variable();
    
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS


    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_points,dim);
    Matrix InvJ0(dim,dim);
    Vector temp(number_of_points);

    Vector heat_flux_local(number_of_points);
    Vector nodal_conductivity_m1(number_of_points);
    Vector nodal_conductivity_m2(number_of_points);
    Vector nodal_conductivity_m3(number_of_points);
    Vector nodal_conductivity_m4(number_of_points);
    Vector nodal_conductivity_m5(number_of_points);
    Vector nodal_conductivity_m6(number_of_points);
    Vector nodal_conductivity_m7(number_of_points);
    Vector nodal_density_m1(number_of_points);
    Vector nodal_density_m2(number_of_points);
    Vector nodal_density_m3(number_of_points);
    Vector nodal_density_m4(number_of_points);
    Vector nodal_density_m5(number_of_points);
    Vector nodal_density_m6(number_of_points);
    Vector nodal_density_m7(number_of_points);
    Vector nodal_specific_heat_m1(number_of_points);
    Vector nodal_specific_heat_m2(number_of_points);
    Vector nodal_specific_heat_m3(number_of_points);
    Vector nodal_specific_heat_m4(number_of_points);
    Vector nodal_specific_heat_m5(number_of_points);
    Vector nodal_specific_heat_m6(number_of_points);
    Vector nodal_specific_heat_m7(number_of_points);
    Vector nodal_mass_fraction_m1(number_of_points);
    Vector nodal_mass_fraction_m2(number_of_points);
    Vector nodal_mass_fraction_m3(number_of_points);
    Vector nodal_mass_fraction_m4(number_of_points);
    Vector nodal_mass_fraction_m5(number_of_points);
    Vector nodal_mass_fraction_m6(number_of_points);
    Vector nodal_mass_fraction_m7(number_of_points);
    
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        heat_flux_local[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_volume_source_var);
        
        nodal_conductivity_m1[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m1_var);
        nodal_conductivity_m2[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m2_var);
        nodal_conductivity_m3[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m3_var);
        nodal_conductivity_m4[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m4_var);
        nodal_conductivity_m5[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m5_var);
        nodal_conductivity_m6[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m6_var);
        nodal_conductivity_m7[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_m7_var);
        nodal_density_m1[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m1_var);
        nodal_density_m2[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m2_var);
        nodal_density_m3[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m3_var);
        nodal_density_m4[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m4_var);
        nodal_density_m5[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m5_var);
        nodal_density_m6[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m6_var);
        nodal_density_m7[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_density_m7_var);
        nodal_specific_heat_m1[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m1_var);
        nodal_specific_heat_m2[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m2_var);
        nodal_specific_heat_m3[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m3_var);
        nodal_specific_heat_m4[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m4_var);
        nodal_specific_heat_m5[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m5_var);
        nodal_specific_heat_m6[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m6_var);
        nodal_specific_heat_m7[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_specific_heat_m7_var);
        nodal_mass_fraction_m1[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m1_var);
        nodal_mass_fraction_m2[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m2_var);
        nodal_mass_fraction_m3[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m3_var);
        nodal_mass_fraction_m4[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m4_var);
        nodal_mass_fraction_m5[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m5_var);
        nodal_mass_fraction_m6[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m6_var);
        nodal_mass_fraction_m7[node_element]  = r_geometry[node_element].FastGetSolutionStepValue(r_mass_fraction_m7_var);
    }

    r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    double DetJ0;
    const double delta_time = r_process_info[DELTA_TIME];
    
    Vector nodal_conductivity_sum(number_of_points);
    Vector nodal_density_sum(number_of_points);
    Vector nodal_specific_heat_sum(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        nodal_conductivity_sum[node_element] = nodal_conductivity_m1[node_element] * nodal_mass_fraction_m1[node_element] + 
                                                                                         nodal_conductivity_m2[node_element] * nodal_mass_fraction_m2[node_element] +
                                                                                         nodal_conductivity_m3[node_element] * nodal_mass_fraction_m3[node_element] +
                                                                                         nodal_conductivity_m4[node_element] * nodal_mass_fraction_m4[node_element] +
                                                                                         nodal_conductivity_m5[node_element] * nodal_mass_fraction_m5[node_element] +
                                                                                         nodal_conductivity_m6[node_element] * nodal_mass_fraction_m6[node_element] +
                                                                                         nodal_conductivity_m7[node_element] * nodal_mass_fraction_m7[node_element] ;
         nodal_density_sum[node_element] = nodal_density_m1[node_element] * nodal_mass_fraction_m1[node_element] + 
                                                                                nodal_density_m2[node_element] * nodal_mass_fraction_m2[node_element] +
                                                                                nodal_density_m3[node_element] * nodal_mass_fraction_m3[node_element] +
                                                                                nodal_density_m4[node_element] * nodal_mass_fraction_m4[node_element] +
                                                                                nodal_density_m5[node_element] * nodal_mass_fraction_m5[node_element] +
                                                                                nodal_density_m6[node_element] * nodal_mass_fraction_m6[node_element] +
                                                                                nodal_density_m7[node_element] * nodal_mass_fraction_m7[node_element] ;
         nodal_specific_heat_sum[node_element] = nodal_specific_heat_m1[node_element] * nodal_mass_fraction_m1[node_element] + 
                                                                                            nodal_specific_heat_m2[node_element] * nodal_mass_fraction_m2[node_element] +
                                                                                            nodal_specific_heat_m3[node_element] * nodal_mass_fraction_m3[node_element] +
                                                                                            nodal_specific_heat_m4[node_element] * nodal_mass_fraction_m4[node_element] +
                                                                                            nodal_specific_heat_m5[node_element] * nodal_mass_fraction_m5[node_element] +
                                                                                            nodal_specific_heat_m6[node_element] * nodal_mass_fraction_m6[node_element] +
                                                                                            nodal_specific_heat_m7[node_element] * nodal_mass_fraction_m7[node_element] ;
    }//Vector有没有重载*运算符
    
    Vector nodal_heat_capacity_sum(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        nodal_heat_capacity_sum[node_element] = nodal_density_sum[node_element] * nodal_specific_heat_sum[node_element]; // ρc
    }
    
    Vector temp_prev_step(number_of_points);
    for (unsigned int i = 0; i < number_of_points; i++)
        temp_prev_step[i] = r_geometry[i].FastGetSolutionStepValue(r_unknown_var, 1); // Previous step temperature

    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point],InvJ0,DetJ0);

        //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point); //these are the N which correspond to the gauss point "i_point"
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double conductivity_gauss = inner_prod(N, nodal_conductivity_sum);
        const double heat_capacity_gauss = inner_prod(N, nodal_heat_capacity_sum); // ρc at gauss point
        
        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); //
        noalias(rLeftHandSideMatrix) += (heat_capacity_gauss / delta_time) * IntToReferenceWeight * outer_prod(N, N); // Mass matrix contribution

        // Calculating the local RHS
        const double qgauss = inner_prod(N, heat_flux_local);

        noalias(rRightHandSideVector) += IntToReferenceWeight*qgauss*N;
        noalias(rRightHandSideVector) += prod((heat_capacity_gauss / delta_time) * IntToReferenceWeight * outer_prod(N, N), temp_prev_step);
    }


    // RHS = ExtForces - K*temp;
    for (unsigned int i = 0; i < number_of_points; i++)
        temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);

    //axpy_prod(rLeftHandSideMatrix, temp, rRightHandSideVector, false);  //RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }
}

//************************************************************************************
//************************************************************************************
int LaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDensityVariable()) << "No Density Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedSpecificHeatVariable()) << "No SpecificHeat Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();
    const Variable<double>& r_density_var = r_settings.GetDensityVariable();
    const Variable<double>& r_specific_heat_var = r_settings.GetSpecificHeatVariable();

    const auto& r_geom = GetGeometry();

    for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
    {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusivity_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_volume_source_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_density_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_specific_heat_var, r_node);

        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
    }

    return Element::Check(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
Element::IntegrationMethod LaplacianElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

} // Namespace Kratos
