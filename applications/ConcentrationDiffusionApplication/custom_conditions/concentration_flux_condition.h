//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//


#if !defined(KRATOS_CONCENTATION_FLUX_CONDITION_H_INCLUDED )
#define  KRATOS_CONCENTRATION_FLUX_CONDITION_H_INCLUDED

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/condition.h"
#include "concentration_diffusion_application_variables.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/variables.h"


namespace Kratos
{

///@addtogroup ConvectionDiffusionApplication
///@{

namespace ConcentrationFluxConditionInternals
{
///@name Auxiliary data structure to hold FEM data
///@{

template< unsigned int TNodeNumber >
class IntegrationData
{
public:

    IntegrationData(
        Geometry< Node >& rGeometry,
        const Variable<double>& rFluxVar
        ):
        mGaussPoint(0),
        mNodalFluxes(TNodeNumber,0.0)
    {
        NumGauss = rGeometry.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_2);
        Vector DetJ = ZeroVector(NumGauss);
        rGeometry.DeterminantOfJacobian(DetJ,GeometryData::IntegrationMethod::GI_GAUSS_2);

        mShapeFunctionValues.resize(NumGauss,TNodeNumber);
        mIntegrationWeights.resize(NumGauss);

        noalias(mShapeFunctionValues) = rGeometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        const auto& IntegrationPoints = rGeometry.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            mIntegrationWeights[g] = DetJ[g]*IntegrationPoints[g].Weight();
        }

        for (unsigned int i = 0; i < TNodeNumber; i++)
        {
            mNodalFluxes[i] = rGeometry[i].FastGetSolutionStepValue(rFluxVar);
        }
    }

    void SetCurrentGaussPoint(unsigned int g)
    {
        mGaussPoint = g;
    }

    double N(unsigned int i) const
    {
        return mShapeFunctionValues(mGaussPoint,i);
    }

    double NodalFlux(unsigned int i) const
    {
        return mNodalFluxes[i];
    }

    double GaussPointFlux() const
    {
        double flux = mNodalFluxes[0]*mShapeFunctionValues(mGaussPoint,0);
        for (unsigned int i = 1; i < TNodeNumber; i++)
        {
            flux += mNodalFluxes[i]*mShapeFunctionValues(mGaussPoint,i);
        }
        return flux;
    }

    double IntegrationWeight() const
    {
        return mIntegrationWeights[mGaussPoint];
    }

    unsigned int NumGauss;

private:

    unsigned int mGaussPoint;

    array_1d<double,TNodeNumber> mNodalFluxes;

    Matrix mShapeFunctionValues;

    Vector mIntegrationWeights;
};

///@}

}

///@}
///@name Kratos Classes
///@{

/// A basic Neumann condition for concentration-diffusion problems.
/** It applies a flux condition based on the nodal values of the
 *  variable defined as SurfaceSourceVariable by the
 *  CONCENTRATION_DIFFUSION_SETTINGS variable in the given ProcessInfo.
 */
template< unsigned int TNodeNumber >
class ConcentrationFluxCondition: public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FluxCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConcentrationFluxCondition);

    typedef Condition::MatrixType MatrixType;
    typedef Condition::VectorType VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    ConcentrationFluxCondition(
        IndexType NewId,
        Geometry< Node >::Pointer pGeometry);

    ConcentrationFluxCondition(
        IndexType NewId,
        Geometry< Node >::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    ~ConcentrationFluxCondition() override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& ConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6 > >& rVariable,
        std::vector<array_1d<double, 6 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:

    ///@name Protected Life cycle
    ///@{

    // Default constructor necessary for serialization
    ConcentrationFluxCondition();

    ///@}
    ///@name Protected Operations
    ///@{

    void AddIntegrationPointRHSContribution(
        VectorType& rRightHandSideVector,
        const ConcentrationFluxConditionInternals::IntegrationData<TNodeNumber>& rData);

    void CalculateNormal(array_1d<double,3>& An );

    ///@}

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ConcentrationFluxCondition& operator=(ConcentrationFluxCondition const& rOther);

    /// Copy constructor.
    ConcentrationFluxCondition(ConcentrationFluxCondition const& rOther);

    ///@}

}; // Class ConcentrationFluxCondition

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONCENTRATION_FLUX_CONDITION_H_INCLUDED  defined
