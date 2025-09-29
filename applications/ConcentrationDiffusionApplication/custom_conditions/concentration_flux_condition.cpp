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

// System includes


// External includes


// Include Base h
#include "custom_conditions/concentration_flux_condition.h"
#include "includes/variables.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
ConcentrationFluxCondition<TNodeNumber>::ConcentrationFluxCondition(IndexType NewId, Geometry< Node >::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

template< unsigned int TNodeNumber >
ConcentrationFluxCondition<TNodeNumber>::ConcentrationFluxCondition(
    IndexType NewId,
    Geometry< Node >::Pointer pGeometry,
    Properties::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

template< unsigned int TNodeNumber >
ConcentrationFluxCondition<TNodeNumber>::~ConcentrationFluxCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
Condition::Pointer ConcentrationFluxCondition<TNodeNumber>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcentrationFluxCondition<TNodeNumber>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TNodeNumber >
Condition::Pointer ConcentrationFluxCondition<TNodeNumber>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcentrationFluxCondition<TNodeNumber>>(NewId, pGeom, pProperties);
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNodeNumber)
    {
        rLeftHandSideMatrix.resize(TNodeNumber,TNodeNumber,false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNodeNumber,TNodeNumber);

    this->CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNodeNumber)
    {
        rRightHandSideVector.resize(TNodeNumber,false);
    }
    noalias(rRightHandSideVector) = ZeroVector(TNodeNumber);

    ConcentrationFluxConditionInternals::IntegrationData<TNodeNumber> Data(this->GetGeometry(), FACE_CONCENTRATION_FLUX);

    for (unsigned int g = 0; g < Data.NumGauss; g++)
    {
        Data.SetCurrentGaussPoint(g);

        this->AddIntegrationPointRHSContribution(rRightHandSideVector, Data);
    }

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != TNodeNumber)
    {
        rResult.resize(TNodeNumber,false);
    }

    const Geometry< Node >& rGeometry = this->GetGeometry();

    for (unsigned int i = 0; i < TNodeNumber; i++)
    {
        rResult[i] = rGeometry[i].GetDof(CONCENTRATION).EquationId();
    }

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rConditionalDofList.size() != TNodeNumber)
    {
        rConditionalDofList.resize(TNodeNumber);
    }

    const Geometry< Node >& rGeometry = this->GetGeometry();

    for (unsigned int i = 0; i < TNodeNumber; i++)
    {
        rConditionalDofList[i] = rGeometry[i].pGetDof(CONCENTRATION);
    }

    KRATOS_CATCH("")
}

template< unsigned int TNodeNumber >
GeometryData::IntegrationMethod ConcentrationFluxCondition<TNodeNumber>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> > &rVariable,
    std::vector<array_1d<double,3> > &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    if (rVariable == NORMAL)
    {
        this->CalculateNormal(rValues[0]);
    }
    else
    {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const ConcentrationFluxCondition* const_this = static_cast< const ConcentrationFluxCondition* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    // Copy the values to the different gauss points
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ConcentrationFluxCondition* const_this = static_cast< const ConcentrationFluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        rValues[g] = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ConcentrationFluxCondition* const_this = static_cast< const ConcentrationFluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ConcentrationFluxCondition* const_this = static_cast< const ConcentrationFluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ConcentrationFluxCondition* const_this = static_cast< const ConcentrationFluxCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

// Input and Output ///////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
std::string ConcentrationFluxCondition<TNodeNumber>::Info() const
{
    std::stringstream buffer;
    buffer << "ConcentrationFluxCondition #" << Id();
    return buffer.str();
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ConcentrationFluxCondition #" << Id();
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::PrintData(std::ostream& rOStream) const
{
    rOStream << "ConcentrationFluxCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}


// Finite element functions ///////////////////////////////////////////////////

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::AddIntegrationPointRHSContribution(
    VectorType& rRightHandSideVector,
    const ConcentrationFluxConditionInternals::IntegrationData<TNodeNumber>& rData)
{
    double InterpolatedFlux = rData.GaussPointFlux();
    for (unsigned int i = 0; i < TNodeNumber; i++)
    {
        rRightHandSideVector[i] += rData.N(i) * InterpolatedFlux * rData.IntegrationWeight();
    }
}


template <>
void ConcentrationFluxCondition<2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;
}


template <>
void ConcentrationFluxCondition<3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}

template <>
void ConcentrationFluxCondition<4>::CalculateNormal(array_1d<double,3>& An )
{
    KRATOS_ERROR << "This function is not yet implemented" << std::endl;

}

// Serialization //////////////////////////////////////////////////////////////

template< unsigned int TNodeNumber >
ConcentrationFluxCondition<TNodeNumber>::ConcentrationFluxCondition():
    Condition()
{
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template< unsigned int TNodeNumber >
void ConcentrationFluxCondition<TNodeNumber>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

template class ConcentrationFluxCondition<2>;
template class ConcentrationFluxCondition<3>;
template class ConcentrationFluxCondition<4>;

}
