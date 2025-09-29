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


// Application includes
#include "axisymmetric_concentration_face.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

AxisymmetricConcentrationFace::AxisymmetricConcentrationFace(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : ConcentrationFace(NewId, pGeometry)
{
}

AxisymmetricConcentrationFace::AxisymmetricConcentrationFace(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    typename PropertiesType::Pointer pProperties)
    : ConcentrationFace(NewId, pGeometry, pProperties)
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer AxisymmetricConcentrationFace::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AxisymmetricConcentrationFace>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer AxisymmetricConcentrationFace::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AxisymmetricConcentrationFace>(NewId, pGeom, pProperties);
}

// Protected Operations //////////////////////////////////////////////////////////

void AxisymmetricConcentrationFace::SetIntegrationWeight(
    const IndexType IntegrationPointIndex,
    const typename GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const Vector &rJacobianDeterminantsVector,
    ConditionDataStruct &rData)
{
    double radius = 0.0;
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < r_geom.PointsNumber(); ++i) {
        radius += rData.N[i] * r_geom[i].Y();
    }
    rData.Weight = 2.0 * Globals::Pi * radius * rJacobianDeterminantsVector[IntegrationPointIndex] * rIntegrationPoints[IntegrationPointIndex].Weight();
}

// Input and Output ///////////////////////////////////////////////////////////

std::string AxisymmetricConcentrationFace::Info() const
{
    std::stringstream buffer;
    buffer << "AxisymmetricConcentrationFace #" << Id();
    return buffer.str();
}

void AxisymmetricConcentrationFace::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "AxisymmetricConcentrationFace #" << Id();
}

void AxisymmetricConcentrationFace::PrintData(std::ostream& rOStream) const
{
    rOStream << "AxisymmetricConcentrationFace #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

AxisymmetricConcentrationFace::AxisymmetricConcentrationFace():
    BaseType()
{
}

void AxisymmetricConcentrationFace::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

void AxisymmetricConcentrationFace::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

}
