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

//~ axis y is symmetry axis

#pragma once

// System includes


// External includes


// Project includes
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/serializer.h"

// Application includes
#include "eulerian_concen_diff.h"

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(CONCENTRATION_DIFFUSION_APPLICATION) AxisymmetricEulerianConcentrationDiffusionElement
    : public EulerianConcentrationDiffusionElement<TDim, TNumNodes>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of AxisymmetricEulerianConcentrationDiffusionElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AxisymmetricEulerianConcentrationDiffusionElement);

    /// Base concentration-diffusion element type
    using BaseType = EulerianConcentrationDiffusionElement<TDim, TNumNodes>;

    /// Geometry type
    using GeometryType = typename BaseType::GeometryType;

    /// Properties type
    using PropertiesType = typename BaseType::PropertiesType;

    /// Index type
    using IndexType = typename BaseType::IndexType;

    /// Size type
    using SizeType = typename BaseType::SizeType;

    /// Vector type
    using VectorType = typename BaseType::VectorType;

    /// Matrix type
    using MatrixType = typename BaseType::MatrixType;

    /// Nodes array type
    using NodesArrayType = typename BaseType::NodesArrayType;

    /// Shape functions gradient container type
    using ShapeFunctionsGradientsType = typename GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructors.
    AxisymmetricEulerianConcentrationDiffusionElement() : BaseType()
    {
    }

    AxisymmetricEulerianConcentrationDiffusionElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {}

    AxisymmetricEulerianConcentrationDiffusionElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~AxisymmetricEulerianConcentrationDiffusionElement(){};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AxisymmetricEulerianConcentrationDiffusionElement<TDim, TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AxisymmetricEulerianConcentrationDiffusionElement<TDim, TNumNodes>>(NewId, pGeom, pProperties);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "AxisymmetricConcentrationDiffusion #";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << this->Id();
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{

    /// Integration rule to be employed
    static constexpr GeometryData::IntegrationMethod mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class AxisymmetricConcentrationDiffusion

///@}

///@} // ConcentrationDiffusionApplication group

} // namespace Kratos.
