// KRATOS MULTIPHYSICS COUPLING ELEMENT
// 
// Author: Your Name
// 
// Thermo-Chemo-Mechanical-PhaseChange Coupled Element Implementation

#include "custom_elements/thermo_diff_struct_phase_element.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/structural_analysis_element_utilities.h"

namespace Kratos
{
/// Default constructor
template <unsigned int TDim, unsigned int TNumNodes>
ThermoDiffStructPhaseElement<TDim, TNumNodes>::ThermoDiffStructPhaseElement() : BaseType() {}

/// Constructor with Id and geometry
template <unsigned int TDim, unsigned int TNumNodes>
ThermoDiffStructPhaseElement<TDim, TNumNodes>::ThermoDiffStructPhaseElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry) {}

/// Constructor with Id, geometry and properties
template <unsigned int TDim, unsigned int TNumNodes>
ThermoDiffStructPhaseElement<TDim, TNumNodes>::ThermoDiffStructPhaseElement(IndexType NewId, 
                                        GeometryType::Pointer pGeometry, 
                                        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties) {}

// copy constructor
template <unsigned int TDim, unsigned int TNumNodes>
ThermoDiffStructPhaseElement<TDim, TNumNodes>::ThermoDiffStructPhaseElement(ThermoDiffStructPhaseElement const& rOther)
    : BaseType(rOther) 
{}

/***********************************************************************************/
/***********************************************************************************/
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer ThermoDiffStructPhaseElement<TDim, TNumNodes>::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    ThermoDiffStructPhaseElement::Pointer p_new_elem = Kratos::make_intrusive<ThermoDiffStructPhaseElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/


template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, 
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
        
        const auto& r_geom = GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        
        // 计算总自由度数量
        SizeType total_dof;
        if (dimension == 2) {
            total_dof = number_of_nodes * DOF_PER_NODE_2D;
        } else {
            total_dof = number_of_nodes * DOF_PER_NODE_3D;
        }
        
        if (rResult.size() != total_dof)
            rResult.resize(total_dof, false);
        
        // 获取各个自由度的位置
        const SizeType disp_x_pos = r_geom[0].GetDofPosition(DISPLACEMENT_X);
        const SizeType temp_pos = r_geom[0].GetDofPosition(TEMPERATURE);
        const SizeType conc_pos = r_geom[0].GetDofPosition(CONCENTRATION);
        
        if (dimension == 2) {
            // 2D情况：每个节点4个自由度 [ux, uy, T, C]
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const SizeType base_index = i * DofsPerNode2D;
                
                // 位移自由度
                rResult[base_index]     = r_geom[i].GetDof(DISPLACEMENT_X, disp_x_pos).EquationId();
                rResult[base_index + 1] = r_geom[i].GetDof(DISPLACEMENT_Y, disp_x_pos + 1).EquationId();
                
                // 温度自由度
                rResult[base_index + 2] = r_geom[i].GetDof(TEMPERATURE, temp_pos).EquationId();
                
                // 浓度自由度
                rResult[base_index + 3] = r_geom[i].GetDof(CONCENTRATION, conc_pos).EquationId();
            }
        } else {
            // 3D情况：每个节点5个自由度 [ux, uy, uz, T, C]
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const SizeType base_index = i * DofsPerNode3D;
                
                // 位移自由度
                rResult[base_index]     = r_geom[i].GetDof(DISPLACEMENT_X, disp_x_pos).EquationId();
                rResult[base_index + 1] = r_geom[i].GetDof(DISPLACEMENT_Y, disp_x_pos + 1).EquationId();
                rResult[base_index + 2] = r_geom[i].GetDof(DISPLACEMENT_Z, disp_x_pos + 2).EquationId();
                
                // 温度自由度
                rResult[base_index + 3] = r_geom[i].GetDof(TEMPERATURE, temp_pos).EquationId();
                
                // 浓度自由度
                rResult[base_index + 4] = r_geom[i].GetDof(CONCENTRATION, conc_pos).EquationId();
            }
        }
        
        KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::GetDofList(
    DofsVectorType& ElementalDofList, 
    const ProcessInfo& rCurrentProcessInfo) const
{
     KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    
    // 计算总自由度数量并直接调整大小
    SizeType total_dof;
    if (dimension == 2) {
        total_dof = number_of_nodes * DofsPerNode2D;
    } else {
        total_dof = number_of_nodes * DofsPerNode3D;
    }
    
    if (rElementalDofList.size() != total_dof)
        rElementalDofList.resize(total_dof);

    if (dimension == 2) {
        // 2D情况：每个节点按 [ux, uy, T, C] 顺序填充
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType base_index = i * DofsPerNode2D;
            
            rElementalDofList[base_index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[base_index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[base_index + 2] = r_geom[i].pGetDof(TEMPERATURE);
            rElementalDofList[base_index + 3] = r_geom[i].pGetDof(CONCENTRATION);
        }
    } else {
        // 3D情况：每个节点按 [ux, uy, uz, T, C] 顺序填充
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType base_index = i * DofsPerNode3D;
            
            rElementalDofList[base_index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[base_index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[base_index + 2] = r_geom[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[base_index + 3] = r_geom[i].pGetDof(TEMPERATURE);
            rElementalDofList[base_index + 4] = r_geom[i].pGetDof(CONCENTRATION);
        }
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize system matrices
    if (rLeftHandSideMatrix.size1() != TotalDofs || rLeftHandSideMatrix.size2() != TotalDofs)
        rLeftHandSideMatrix.resize(TotalDofs, TotalDofs, false);
    
    if (rRightHandSideVector.size() != TotalDofs)
        rRightHandSideVector.resize(TotalDofs, false);

    // Initialize to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(TotalDofs, TotalDofs);
    noalias(rRightHandSideVector) = ZeroVector(TotalDofs);

    // Initialize element variables
    CoupledElementVariables Variables;
    InitializeCoupledElement(Variables, rCurrentProcessInfo);
    
    // Get nodal values
    GetAllNodalValues(Variables, rCurrentProcessInfo);

    // Integration loop
    const auto& integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        // Calculate geometry data
        CalculateGeometryData(Variables, point_number);
        
        // Update phase transformations for each node
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
        {
            auto phase_results = CalculatePhaseTransformation(Variables, i_node, rCurrentProcessInfo);
            
            // Update latent heat and transformation strain
            Variables.latent_heat_generation[i_node] = phase_results.total_latent_heat;
            
            // Store phase transformation strain (simplified - use average)
            if (phase_results.transformation_occurred)
            {
                for (IndexType i = 0; i < 6; ++i)
                {
                    Variables.phase_transformation_strain[i][i_node] = phase_results.phase_strain_increment[i];
                }
            }
        }
        
        // Update material properties based on current phase composition
        UpdateMaterialProperties(Variables, rCurrentProcessInfo);
        
        // Calculate individual field contributions
        CalculateThermalContribution(rLeftHandSideMatrix, rRightHandSideVector, Variables, rCurrentProcessInfo);
        CalculateDiffusionContribution(rLeftHandSideMatrix, rRightHandSideVector, Variables, rCurrentProcessInfo);
        CalculateMechanicalContribution(rLeftHandSideMatrix, rRightHandSideVector, Variables, rCurrentProcessInfo);
        
        // Add coupling terms
        AddCouplingTerms(rLeftHandSideMatrix, rRightHandSideVector, Variables, rCurrentProcessInfo);
    }

    KRATOS_CATCH("Error in ThermoDiffStructPhaseElement::CalculateLocalSystem")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp_matrix;
    CalculateLocalSystem(temp_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::InitializeCoupledElement(
    CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    rVariables.dt = rCurrentProcessInfo[DELTA_TIME];
    rVariables.dt_inv = 1.0 / rVariables.dt;
    rVariables.lumping_factor = 1.0 / double(TNumNodes);

    // Initialize material properties
    rVariables.effective_density = 0.0;
    rVariables.effective_specific_heat = 0.0;
    rVariables.effective_conductivity = 0.0;
    rVariables.effective_diffusivity = 0.0;
    rVariables.effective_young_modulus = 0.0;
    rVariables.effective_poisson_ratio = 0.0;

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::GetAllNodalValues(
    CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto& r_properties = GetProperties();
    
    bool is_multiphase_defined = false;
    if (r_properties.Has(IS_DEFINED_MULTIPHASE))
        is_multiphase_defined = r_properties[IS_DEFINED_MULTIPHASE];

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        // Temperature values
        rVariables.temperature_current[i] = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE, 0);
        rVariables.temperature_previous[i] = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE, 1);
        
        // Concentration values
        rVariables.concentration_current[i] = r_geometry[i].FastGetSolutionStepValue(CONCENTRATION, 0);
        rVariables.concentration_previous[i] = r_geometry[i].FastGetSolutionStepValue(CONCENTRATION, 1);
        
        // Displacement values
        for (IndexType j = 0; j < TDim; ++j)
        {
            rVariables.displacement_current[j][i] = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_X + j, 0);
            rVariables.displacement_previous[j][i] = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_X + j, 1);
        }

        // Phase fractions
        if (is_multiphase_defined)
        {
            for (IndexType phase = 0; phase < NumPhases; ++phase)
            {
                // Assuming variables MASS_FRACTION_M1 through MASS_FRACTION_M7
                rVariables.phase_fractions_current[i][phase] = 
                    r_geometry[i].FastGetSolutionStepValue(MASS_FRACTION_M1 + phase, 0);
                rVariables.phase_fractions_previous[i][phase] = 
                    r_geometry[i].FastGetSolutionStepValue(MASS_FRACTION_M1 + phase, 1);
            }
        }
        else
        {
            // Initialize to single phase
            rVariables.phase_fractions_current[i].fill(0.0);
            rVariables.phase_fractions_previous[i].fill(0.0);
            rVariables.phase_fractions_current[i][0] = 1.0;  // Assume first phase is 100%
            rVariables.phase_fractions_previous[i][0] = 1.0;
        }

        // Source terms
        rVariables.thermal_source[i] = r_geometry[i].FastGetSolutionStepValue(HEAT_FLUX, 0);
        // Assuming you have a concentration source variable
        // rVariables.concentration_source[i] = r_geometry[i].FastGetSolutionStepValue(CONCENTRATION_SOURCE, 0);
        
        // Body forces
        for (IndexType j = 0; j < TDim; ++j)
        {
            rVariables.body_force[j][i] = r_geometry[i].FastGetSolutionStepValue(BODY_FORCE_X + j, 0);
        }
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateGeometryData(
    CoupledElementVariables& rVariables,
    const IndexType IntegrationPoint)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto& integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    
    // Get shape functions
    rVariables.N = r_geometry.ShapeFunctionsValues(rVariables.N, integration_points[IntegrationPoint].Coordinates());
    
    // Calculate derivatives
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0, GetIntegrationMethod());
    
    Matrix InvJ0;
    double detJ0;
    MathUtils<double>::InvertMatrix(J0[IntegrationPoint], InvJ0, detJ0);
    
    rVariables.detJ = detJ0;
    
    GeometryType::ShapeFunctionsGradientsType DN_De;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_De, GetIntegrationMethod());
    
    noalias(rVariables.DN_DX) = prod(DN_De[IntegrationPoint], InvJ0);
    
    // Integration weight
    rVariables.integration_weight = integration_points[IntegrationPoint].Weight() * detJ0;
    
    // Apply thickness for 2D problems
    if constexpr (TDim == 2)
    {
        if (GetProperties().Has(THICKNESS))
            rVariables.integration_weight *= GetProperties()[THICKNESS];
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::UpdateMaterialProperties(
    CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    const auto& r_geometry = GetGeometry();
    
    rVariables.effective_density = 0.0;
    rVariables.effective_specific_heat = 0.0;
    rVariables.effective_conductivity = 0.0;
    rVariables.effective_diffusivity = 0.0;
    rVariables.effective_young_modulus = 0.0;
    rVariables.effective_poisson_ratio = 0.0;

    // Average over nodes using shape functions
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        double node_density = 0.0;
        double node_specific_heat = 0.0;
        double node_conductivity = 0.0;
        double node_diffusivity = 0.0;
        double node_young_modulus = 0.0;
        double node_poisson_ratio = 0.0;

        // Mix properties based on phase fractions
        for (IndexType phase = 0; phase < NumPhases; ++phase)
        {
            double phase_fraction = rVariables.phase_fractions_current[i][phase];
            
            if (phase_fraction > 1e-12)  // Only consider phases with significant fraction
            {
                // Get phase-specific properties (these need to be defined in your application)
                node_density += phase_fraction * r_properties.GetValue(DENSITY_M1 + phase, r_geometry[i], rCurrentProcessInfo);
                node_specific_heat += phase_fraction * r_properties.GetValue(SPECIFIC_HEAT_M1 + phase, r_geometry[i], rCurrentProcessInfo);
                node_conductivity += phase_fraction * r_properties.GetValue(CONDUCTIVITY_M1 + phase, r_geometry[i], rCurrentProcessInfo);
                // Add other properties as needed
            }
        }

        // Weight by shape function
        rVariables.effective_density += rVariables.N[i] * node_density;
        rVariables.effective_specific_heat += rVariables.N[i] * node_specific_heat;
        rVariables.effective_conductivity += rVariables.N[i] * node_conductivity;
        rVariables.effective_diffusivity += rVariables.N[i] * node_diffusivity;
        rVariables.effective_young_modulus += rVariables.N[i] * node_young_modulus;
        rVariables.effective_poisson_ratio += rVariables.N[i] * node_poisson_ratio;
    }

    // Ensure non-zero values
    if (rVariables.effective_density < 1e-12) rVariables.effective_density = 1.0;
    if (rVariables.effective_specific_heat < 1e-12) rVariables.effective_specific_heat = 1.0;
    if (rVariables.effective_conductivity < 1e-12) rVariables.effective_conductivity = 1.0;

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
typename ThermoDiffStructPhaseElement<TDim, TNumNodes>::PhaseTransformationResults
ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculatePhaseTransformation(
    const CoupledElementVariables& rVariables,
    const IndexType NodeIndex,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    PhaseTransformationResults results;
    results.transformation_occurred = false;
    results.total_latent_heat = 0.0;
    results.phase_fraction_increments.fill(0.0);
    results.reaction_rates.fill(0.0);
    
    if constexpr (TDim == 3) {
        results.phase_strain_increment = ZeroVector(6);
    } else {
        results.phase_strain_increment = ZeroVector(3);
    }

    double temperature = rVariables.temperature_current[NodeIndex];
    double concentration = rVariables.concentration_current[NodeIndex];
    const auto& current_phases = rVariables.phase_fractions_current[NodeIndex];

    // Calculate transformation rates for each transformation type
    results.reaction_rates[0] = CalculateAusteniteTransformationRate(temperature, concentration, current_phases);
    results.reaction_rates[1] = CalculateEutectoidTransformationRate(temperature, concentration, current_phases);
    results.reaction_rates[2] = CalculateBainiteTransformationRate(temperature, concentration, current_phases);
    results.reaction_rates[3] = CalculateMartensiteTransformationRate(temperature, concentration, current_phases);

    // Check if any transformation is significant
    double total_rate = 0.0;
    for (const auto& rate : results.reaction_rates) {
        total_rate += std::abs(rate);
    }

    if (total_rate > 1e-12) {
        results.transformation_occurred = true;
        
        // Calculate phase fraction changes based on transformation kinetics
        // This is a simplified example - you need to implement your specific kinetics
        double dt = rVariables.dt;
        
        // Example: Simple rate-based updates
        results.phase_fraction_increments[1] = results.reaction_rates[0] * dt; // Austenite
        results.phase_fraction_increments[2] = results.reaction_rates[1] * dt; // Eutectoid
        results.phase_fraction_increments[3] = results.reaction_rates[2] * dt; // Bainite
        results.phase_fraction_increments[4] = results.reaction_rates[3] * dt; // Martensite
        
        // Calculate latent heat from transformations
        const auto& r_properties = GetProperties();
        const auto& r_geometry = GetGeometry();
        
        results.total_latent_heat += r_properties.GetValue(LATENT_AUSTENIZE) * results.phase_fraction_increments[1];
        results.total_latent_heat += r_properties.GetValue(LATENT_EUTECTOID_DECOMPOSITION) * results.phase_fraction_increments[2];
        results.total_latent_heat += r_properties.GetValue(LATENT_BAINITE_TRANSFORMATION) * results.phase_fraction_increments[3];
        results.total_latent_heat += r_properties.GetValue(LATENT_MARTENSITE_TRANSFORMATION) * results.phase_fraction_increments[4];
        
        // Calculate phase transformation strain (simplified volumetric strain)
        double volumetric_strain = 0.0;
        for (IndexType i = 1; i < 5; ++i) {
            // Assume each phase has different specific volume
            double phase_volumetric_strain = results.phase_fraction_increments[i] * 0.01; // 1% strain per phase change
            volumetric_strain += phase_volumetric_strain;
        }
        
        // Apply volumetric strain to strain tensor
        if constexpr (TDim == 3) {
            results.phase_strain_increment[0] = volumetric_strain / 3.0; // �xx
            results.phase_strain_increment[1] = volumetric_strain / 3.0; // �yy  
            results.phase_strain_increment[2] = volumetric_strain / 3.0; // �zz
            // Shear strains remain zero for volumetric transformation
        } else {
            results.phase_strain_increment[0] = volumetric_strain / 2.0; // �xx
            results.phase_strain_increment[1] = volumetric_strain / 2.0; // �yy
        }
    }

    return results;

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateThermalContribution(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Thermal capacity matrix
    BoundedMatrix<double, TNumNodes, TNumNodes> thermal_mass_matrix = ZeroMatrix(TNumNodes, TNumNodes);
    
    // Conductivity matrix  
    BoundedMatrix<double, TNumNodes, TNumNodes> thermal_stiffness_matrix = ZeroMatrix(TNumNodes, TNumNodes);
    
    // Build thermal mass matrix (capacity matrix)
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            thermal_mass_matrix(i, j) = rVariables.effective_density * rVariables.effective_specific_heat * 
                                       rVariables.N[i] * rVariables.N[j] * rVariables.integration_weight;
        }
    }
    
    // Build thermal stiffness matrix (conductivity matrix)
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            double conductivity_term = 0.0;
            for (IndexType k = 0; k < TDim; ++k) {
                conductivity_term += rVariables.DN_DX(i, k) * rVariables.DN_DX(j, k);
            }
            thermal_stiffness_matrix(i, j) = rVariables.effective_conductivity * conductivity_term * 
                                           rVariables.integration_weight;
        }
    }

    // Assembly into global matrix (temperature DOFs)
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            IndexType global_i = GetDofIndex(i, TDim);     // Temperature DOF
            IndexType global_j = GetDofIndex(j, TDim);     // Temperature DOF
            
            // Transient term
            rLeftHandSideMatrix(global_i, global_j) += rVariables.dt_inv * thermal_mass_matrix(i, j);
            
            // Diffusion term
            rLeftHandSideMatrix(global_i, global_j) += thermal_stiffness_matrix(i, j);
        }
        
        // RHS: Previous time step + sources
        IndexType global_i = GetDofIndex(i, TDim);
        
        // Previous temperature contribution
        for (IndexType j = 0; j < TNumNodes; ++j) {
            rRightHandSideVector[global_i] += rVariables.dt_inv * thermal_mass_matrix(i, j) * 
                                            rVariables.temperature_previous[j];
        }
        
        // Heat sources (including latent heat)
        rRightHandSideVector[global_i] += (rVariables.thermal_source[i] + rVariables.latent_heat_generation[i]) * 
                                        rVariables.N[i] * rVariables.integration_weight;
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateDiffusionContribution(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Similar to thermal, but for concentration
    BoundedMatrix<double, TNumNodes, TNumNodes> diffusion_mass_matrix = ZeroMatrix(TNumNodes, TNumNodes);
    BoundedMatrix<double, TNumNodes, TNumNodes> diffusion_stiffness_matrix = ZeroMatrix(TNumNodes, TNumNodes);
    
    // Build diffusion mass matrix
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            diffusion_mass_matrix(i, j) = rVariables.N[i] * rVariables.N[j] * rVariables.integration_weight;
        }
    }
    
    // Build diffusion stiffness matrix
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            double diffusivity_term = 0.0;
            for (IndexType k = 0; k < TDim; ++k) {
                diffusivity_term += rVariables.DN_DX(i, k) * rVariables.DN_DX(j, k);
            }
            diffusion_stiffness_matrix(i, j) = rVariables.effective_diffusivity * diffusivity_term * 
                                             rVariables.integration_weight;
        }
    }

    // Assembly into global matrix (concentration DOFs)
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            IndexType global_i = GetDofIndex(i, TDim + 1); // Concentration DOF
            IndexType global_j = GetDofIndex(j, TDim + 1); // Concentration DOF
            
            // Transient term
            rLeftHandSideMatrix(global_i, global_j) += rVariables.dt_inv * diffusion_mass_matrix(i, j);
            
            // Diffusion term
            rLeftHandSideMatrix(global_i, global_j) += diffusion_stiffness_matrix(i, j);
        }
        
        // RHS: Previous time step + sources
        IndexType global_i = GetDofIndex(i, TDim + 1);
        
        // Previous concentration contribution
        for (IndexType j = 0; j < TNumNodes; ++j) {
            rRightHandSideVector[global_i] += rVariables.dt_inv * diffusion_mass_matrix(i, j) * 
                                            rVariables.concentration_previous[j];
        }
        
        // Concentration sources
        rRightHandSideVector[global_i] += rVariables.concentration_source[i] * rVariables.N[i] * 
                                        rVariables.integration_weight;
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateMechanicalContribution(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Use B-matrix approach similar to SmallDisplacement element
    const IndexType strain_size = (TDim == 3) ? 6 : 3;
    Matrix B_matrix = ZeroMatrix(strain_size, TNumNodes * TDim);
    
    // Build B matrix (strain-displacement relationship)
    StructuralAnalysisElementUtilities::CalculateB(*this, rVariables.DN_DX, B_matrix);
    
    // Constitutive matrix (simplified elastic)
    Matrix D_matrix = ZeroMatrix(strain_size, strain_size);
    
    // Build elastic constitutive matrix
    double E = rVariables.effective_young_modulus;
    double nu = rVariables.effective_poisson_ratio;
    
    if constexpr (TDim == 3) {
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D_matrix(0, 0) = D_matrix(1, 1) = D_matrix(2, 2) = factor * (1.0 - nu);
        D_matrix(0, 1) = D_matrix(1, 0) = D_matrix(0, 2) = D_matrix(2, 0) = 
        D_matrix(1, 2) = D_matrix(2, 1) = factor * nu;
        D_matrix(3, 3) = D_matrix(4, 4) = D_matrix(5, 5) = factor * (1.0 - 2.0 * nu) / 2.0;
    } else {
        // Plane strain
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D_matrix(0, 0) = D_matrix(1, 1) = factor * (1.0 - nu);
        D_matrix(0, 1) = D_matrix(1, 0) = factor * nu;
        D_matrix(2, 2) = factor * (1.0 - 2.0 * nu) / 2.0;
    }

    // Mechanical stiffness matrix
    Matrix mechanical_stiffness = prod(trans(B_matrix), Matrix(prod(D_matrix, B_matrix))) * 
                                 rVariables.integration_weight;
    
    // Assembly mechanical stiffness into global matrix
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType a = 0; a < TDim; ++a) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                for (IndexType b = 0; b < TDim; ++b) {
                    IndexType global_i = GetDofIndex(i, a);
                    IndexType global_j = GetDofIndex(j, b);
                    IndexType local_i = i * TDim + a;
                    IndexType local_j = j * TDim + b;
                    
                    rLeftHandSideMatrix(global_i, global_j) += mechanical_stiffness(local_i, local_j);
                }
            }
        }
    }

    // Body forces and phase transformation stresses
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType a = 0; a < TDim; ++a) {
            IndexType global_i = GetDofIndex(i, a);
            
            // Body forces
            rRightHandSideVector[global_i] += rVariables.body_force[a][i] * rVariables.N[i] * 
                                            rVariables.integration_weight;
            
            // Phase transformation induced stresses (simplified)
            Vector phase_strain = ZeroVector(strain_size);
            for (IndexType k = 0; k < strain_size; ++k) {
                phase_strain[k] = rVariables.phase_transformation_strain[k][i];
            }
            
            Vector phase_stress = prod(D_matrix, phase_strain);
            Vector internal_forces = prod(trans(B_matrix), phase_stress) * rVariables.integration_weight;
            
            rRightHandSideVector[global_i] -= internal_forces[i * TDim + a];
        }
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void ThermoDiffStructPhaseElement<TDim, TNumNodes>::AddCouplingTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const CoupledElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Thermal expansion coupling (Temperature -> Displacement)
    double thermal_expansion_coeff = 1e-5; // Should come from properties
    
    for (IndexType i = 0; i < TNumNodes; ++i) {
        double temperature_change = rVariables.temperature_current[i] - rVariables.temperature_previous[i];
        
        for (IndexType a = 0; a < TDim; ++a) {
            IndexType disp_dof = GetDofIndex(i, a);
            
            // Thermal expansion strain contribution
            double thermal_strain = thermal_expansion_coeff * temperature_change;
            rRightHandSideVector[disp_dof] += thermal_strain * rVariables.N[i] * rVariables.integration_weight;
        }
    }
    
    // Chemical expansion coupling (Concentration -> Displacement)  
    double chemical_expansion_coeff = 1e-6; // Should come from properties
    
    for (IndexType i = 0; i < TNumNodes; ++i) {
        double concentration_change = rVariables.concentration_current[i] - rVariables.concentration_previous[i];
        
        for (IndexType a = 0; a < TDim; ++a) {
            IndexType disp_dof = GetDofIndex(i, a);
            
            // Chemical expansion strain contribution
            double chemical_strain = chemical_expansion_coeff * concentration_change;
            rRightHandSideVector[disp_dof] += chemical_strain * rVariables.N[i] * rVariables.integration_weight;
        }
    }

    // Deformation heating (Displacement -> Temperature)
    // This requires calculating plastic work or deformation energy

    KRATOS_CATCH("")
}

// Private transformation rate calculation methods
template<unsigned int TDim, unsigned int TNumNodes>
double ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateAusteniteTransformationRate(
    double temperature, double concentration, const std::array<double, NumPhases>& current_phases)
{
    // Implement your specific austenite transformation kinetics
    // This is a placeholder - replace with actual model
    double austenite_start_temp = 900.0; // �C
    double rate = 0.0;
    
    if (temperature > austenite_start_temp && current_phases[1] < 0.99) {
        double driving_force = temperature - austenite_start_temp;
        rate = 0.01 * driving_force * (1.0 - current_phases[1]); // Simple kinetics
    }
    
    return rate;
}

template<unsigned int TDim, unsigned int TNumNodes>
double ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateEutectoidTransformationRate(
    double temperature, double concentration, const std::array<double, NumPhases>& current_phases)
{
    // Implement eutectoid transformation kinetics
    double eutectoid_temp = 727.0; // �C for steel
    double rate = 0.0;
    
    if (temperature < eutectoid_temp && current_phases[2] < 0.99) {
        double driving_force = eutectoid_temp - temperature;
        rate = 0.005 * driving_force * current_phases[1]; // Rate proportional to austenite
    }
    
    return rate;
}

template<unsigned int TDim, unsigned int TNumNodes>
double ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateBainiteTransformationRate(
    double temperature, double concentration, const std::array<double, NumPhases>& current_phases)
{
    // Implement bainite transformation kinetics
    double bainite_start_temp = 550.0; // �C
    double bainite_finish_temp = 400.0; // �C
    double rate = 0.0;
    
    if (temperature < bainite_start_temp && temperature > bainite_finish_temp && current_phases[3] < 0.99) {
        double driving_force = bainite_start_temp - temperature;
        rate = 0.003 * driving_force * current_phases[1];
    }
    
    return rate;
}

template<unsigned int TDim, unsigned int TNumNodes>
double ThermoDiffStructPhaseElement<TDim, TNumNodes>::CalculateMartensiteTransformationRate(
    double temperature, double concentration, const std::array<double, NumPhases>& current_phases)
{
    // Implement martensite transformation kinetics
    double martensite_start_temp = 300.0; // C - depends on composition
    double rate = 0.0;
    
    if (temperature < martensite_start_temp && current_phases[4] < 0.99) {
        double driving_force = martensite_start_temp - temperature;
        rate = 0.1 * driving_force * current_phases[1]; // Fast transformation
    }
    
    return rate;
}

// Template instantiations
template class ThermoDiffStructPhaseElement<2, 3>; // 2D Triangle
template class ThermoDiffStructPhaseElement<2, 4>; // 2D Quadrilateral  
template class ThermoDiffStructPhaseElement<3, 4>; // 3D Tetrahedron
template class ThermoDiffStructPhaseElement<3, 8>; // 3D Hexahedron

} // namespace Kratos