// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/eulerian_concen_diff.h"
#include "concentration_diffusion_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "axisymmetric_eulerian_concentration_diffusion.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
void AxisymmetricEulerianConcentrationDiffusionElement<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize of LHS and RHS arrays
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    // Initialize LHS and RHS arrays
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Initialize element data container
    typename BaseType::ElementVariables Variables;
    this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

    // Fill element data container with nodal data
    this->GetNodalValues(Variables, rCurrentProcessInfo);

    // Calculate kinematics
    Vector det_J_vect;
    ShapeFunctionsGradientsType DN_DX;
    const auto& r_geom = this->GetGeometry();
    const auto N = r_geom.ShapeFunctionsValues(mIntegrationMethod);
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mIntegrationMethod);

    // Gauss points loop
    double y_g;
    array_1d<double,TNumNodes> N_g;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX_g;
    const auto integration_points = r_geom.IntegrationPoints(mIntegrationMethod);
    const SizeType n_gauss = integration_points.size();
    for (IndexType g = 0; g < n_gauss; ++g) {
        // Get Gauss point data
        noalias(N_g) = row(N, g);
        noalias(DN_DX_g) = DN_DX[g];

        // Calculate Gauss point values
        y_g = 0.0;
        for (IndexType i = 0; i < TNumNodes; ++i) {
            // Gauss point radius
            y_g += N_g[i] * r_geom[i].Y();
        }

        // Calculate axisymmetric integration weight
        const double w_g = 2.0 * Globals::Pi * y_g * integration_points[g].Weight() * det_J_vect[g];

        // Assemble Gauss point LHS and RHS contributions
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                // Source term
                const double aux_source = w_g * N_g[i] * N_g[j];
                rRightHandSideVector(i) += aux_source * Variables.theta * Variables.volumetric_source[j];
                rRightHandSideVector(i) += aux_source * (1.0 - Variables.theta) * Variables.volumetric_source[j]; // In case we make the body force time dependent

                // Dynamic term
                //const double aux_dyn = w_g * Variables.density * Variables.specific_heat * Variables.dt_inv * N_g[i] * N_g[j];
                //密度*比热
                const double aux_dyn = w_g * Variables.dt_inv * N_g[i] * N_g[j];
                rLeftHandSideMatrix(i, j) += aux_dyn;
                rRightHandSideVector(i) -= aux_dyn * (Variables.c_now[j] - Variables.c_old[j]);

                // Diffusive term
                const double aux_diff = w_g * Variables.diffuse_conductivity * (DN_DX_g(i, 0) * DN_DX_g(j, 0) + DN_DX_g(i, 1) * DN_DX_g(j, 1));
                rLeftHandSideMatrix(i, j) += aux_diff * Variables.theta;
                rRightHandSideVector(i) -= aux_diff * Variables.theta * Variables.c_now[j];
                rRightHandSideVector(i) -= aux_diff * (1.0 - Variables.theta) * Variables.c_old[j];
            }
        }
    }

    KRATOS_CATCH("")
}

template< unsigned int TDim, unsigned int TNumNodes >
void AxisymmetricEulerianConcentrationDiffusionElement< TDim, TNumNodes >::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize of RHS array
    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    // Initialize RHS array
    rRightHandSideVector.clear();

    // Initialize element data container
    typename BaseType::ElementVariables Variables;
    this->InitializeEulerianElement(Variables,rCurrentProcessInfo);

    // Fill element data container with nodal data
    this->GetNodalValues(Variables, rCurrentProcessInfo);

    // Calculate kinematics
    Vector det_J_vect;
    ShapeFunctionsGradientsType DN_DX;
    const auto& r_geom = this->GetGeometry();
    const auto N = r_geom.ShapeFunctionsValues(mIntegrationMethod);
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, det_J_vect, mIntegrationMethod);

    // Gauss points loop
    double y_g;
    array_1d<double,TNumNodes> N_g;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX_g;
    const auto integration_points = r_geom.IntegrationPoints(mIntegrationMethod);
    const SizeType n_gauss = integration_points.size();
    for (IndexType g = 0; g < n_gauss; ++g) {
        // Get Gauss point data
        noalias(N_g) = row(N, g);
        noalias(DN_DX_g) = DN_DX[g];

        // Calculate Gauss point values
        y_g = 0.0;
        for (IndexType i = 0; i < TNumNodes; ++i) {
            // Gauss point radius
            y_g += N_g[i] * r_geom[i].Y();
        }

        // Calculate axisymmetric integration weight
        const double w_g = 2.0 * Globals::Pi * y_g * integration_points[g].Weight() * det_J_vect[g];

        // Assemble Gauss point LHS and RHS contributions
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                // Source term
                const double aux_source = w_g * N_g[i] * N_g[j];
                rRightHandSideVector(i) += aux_source * Variables.theta * Variables.volumetric_source[j];
                rRightHandSideVector(i) += aux_source * (1.0 - Variables.theta) * Variables.volumetric_source[j]; // In case we make the body force time dependent

                // Dynamic term
                const double aux_dyn = w_g * Variables.dt_inv * N_g[i] * N_g[j];
                rRightHandSideVector(i) -= aux_dyn * (Variables.c_now[j] - Variables.c_old[j]);

                // Diffusive terms
                const double aux_diff_1 = w_g * Variables.diffuse_conductivity * N_g[i] * DN_DX_g(j,1) / y_g;
                rRightHandSideVector(i) += aux_diff_1 * Variables.theta * Variables.c_now[j];
                rRightHandSideVector(i) += aux_diff_1 * (1.0 - Variables.theta) * Variables.c_old[j];

                const double aux_diff_2 = w_g * Variables.diffuse_conductivity * (DN_DX_g(i, 0) * DN_DX_g(j, 0) + DN_DX_g(i, 1) * DN_DX_g(j, 1));
                rRightHandSideVector(i) -= aux_diff_2 * Variables.theta * Variables.c_now[j];
                rRightHandSideVector(i) -= aux_diff_2 * (1.0 - Variables.theta) * Variables.c_old[j];
            }
        }
    }
}

template< unsigned int TDim, unsigned int TNumNodes >
int AxisymmetricEulerianConcentrationDiffusionElement< TDim, TNumNodes >::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    // Base element check
    int out = BaseType::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Check that there are no negative y-coordinates (radius is always positive)
    const auto& r_geom = this->GetGeometry();
    for (const auto& r_node : r_geom) {
        KRATOS_ERROR_IF(r_node.Y() < 0.0) << "Negative y-coordinate found in node " << r_node.Id() << ". Axisymmetric radius must be positive." << std::endl;
    }

    return 0;
}

template class AxisymmetricEulerianConcentrationDiffusionElement<2,3>;
template class AxisymmetricEulerianConcentrationDiffusionElement<2,4>;

}
