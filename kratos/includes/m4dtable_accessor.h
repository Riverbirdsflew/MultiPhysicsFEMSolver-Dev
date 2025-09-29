//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marx Xu
//

# pragma once

// System includes

// External includes

// Project includes
#include "includes/accessor.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class M4dTableAccessor
 * @ingroup Kratos Core
 * @brief This class defines the way a certain property is accessed according to a table.
 * @brief The tables are supposed to relate double <-> double type entities. 
 * @brief The input variable is suposed to be a nodally accesible one (either historical or not)
 * @author Marx Xu
 */
class KRATOS_API(KRATOS_CORE) M4dTableAccessor : public Accessor
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    /// Variable type
    using VariableType = Variable<double>;

    /// BaseType
    using BaseType = Accessor;

    using DataLocation = Globals::DataLocation;

    using SizeType = std::size_t;

    using VariableListType = std::array<VariableType*, 4>;

    /// Pointer definition of TableAccessor
    KRATOS_CLASS_POINTER_DEFINITION(M4dTableAccessor);

    ///@}
    ///@name Life Cycle
    ///@{
    M4dTableAccessor() {}

    /// Custom constructor
    M4dTableAccessor(SizeType& rDim, VariableType** rpInputVariable, const std::string& rInputVariableType = "node_historical") 
    {
        // We initialize the variable type only once
        if (rInputVariableType == "node_historical") {
            mInputVariableType = Globals::DataLocation::NodeHistorical;
        } else if (rInputVariableType == "node_non_historical") {
            mInputVariableType = Globals::DataLocation::NodeNonHistorical;
        } else if (rInputVariableType == "element") {
            mInputVariableType = Globals::DataLocation::Element;
        } else {
            KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : 'node_historical', 'node_non_historical' and 'element'" << std::endl;
        }

        // check the dimension
        if (rDim > 0 && rDim <= 4) {
            mDim = rDim;

            for(SizeType i = 0; i < mDim; i++)
                mpInputVariable[i] = rpInputVariable[i];
        } else {
            KRATOS_ERROR << "The table_input_variable_dimension is incorrect. Dimension must be between 1 and 4" << std::endl;
        }
    }

    /// Copy constructor
    M4dTableAccessor(const M4dTableAccessor& rOther) 
    : BaseType(rOther), mDim(rOther.mDim), mInputVariableType(rOther.mInputVariableType)
    {
        for(SizeType i = 0; i < mDim; i++)
            mpInputVariable[i] = rOther.mpInputVariable[i];
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Custom method to retrieve double type properties
     * @param rVariables The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rGeometry The geometry considered
     * @param rShapeFunctionVector The shape function of the GP considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetValue(
        const Variable<double>& rVariables,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief Custom method to retrieve double type properties
     * @param rVariables The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rNode The node considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetValue(
        const Variable<double>& rVariables,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief Custom method to retrieve the derivative of double type properties
     * @param rDerVariables The derivative variable considered (double type properties)
     * @param rVariables The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rGeometry The geometry considered
     * @param rShapeFunctionVector The shape function of the GP considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetDerivative(
        const Variable<double>& rDerVariable,
        const Variable<double>& rVariables,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief Custom method to retrieve the derivative of double type properties
     * @param rDerVariables The derivative variable considered (double type properties)
     * @param rVariables The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rNode The node considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetDerivative(
        const Variable<double>& rDerVariable,
        const Variable<double>& rVariables,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief Custom method to retrieve the integration of double type properties
     * @param rIntVariables The integration variable considered (double type properties)
     * @param rVariables The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rGeometry The geometry considered
     * @param rShapeFunctionVector The shape function of the GP considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetIntegration(
        const Variable<double>& rIntVariable,
        const Variable<double>& rVariables,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief Custom method to retrieve the integration of double type properties
     * @param rIntVariables The integration variable considered (double type properties)
     * @param rVariables The variable considered (double type properties)
     * @param rProperties The properties considered
     * @param rNode The node considered
     * @param rProcessInfo The process info considered
     * @return The double type properties
     */
    double GetIntegration(
        const Variable<double>& rIntVariable,
        const Variable<double>& rVariables,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo
        ) const override;

    /**
     * @brief This computes a material property according to a certain
     * nodal Variable<double> m4dtable
     */
    double GetValueFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief This computes a material property according to a certain
     * nodal Variable<double> m4dtable
     */
    double GetValueFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief This computes a derivative of material property according to a certain
     * nodal Variable<double> m4dtable
     */
    double GetDerivativeFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDerIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief This computes a derivative of material property according to a certain
     * nodal Variable<double> m4dtable
     */
    double GetDerivativeFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDerIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief This computes a integration of material property according to a certain
     * nodal Variable<double> m4dtable
     */
    double GetIntegrationFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rIntIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief This computes a integration of material property according to a certain
     * nodal Variable<double> m4dtable
     */
    double GetIntegrationFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rIntIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief Returns the member input variable
     */
    VariableType& GetInputVariable(SizeType const& m) const
    {
        return (m < 4)? *mpInputVariable[m] : *mpInputVariable[0];
    }

    // Getting a pointer to the class
    Accessor::UniquePointer Clone() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "M4dTableAccessor" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "M4dTableAccessor";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {rOStream << "M4dTableAccessor class";}

    ///@}

private:

    ///@name Member Variables
    ///@{

    SizeType mDim = 0;
    VariableListType mpInputVariable = {nullptr, nullptr, nullptr, nullptr};
    Globals::DataLocation mInputVariableType = Globals::DataLocation::NodeHistorical; // NodalHistorical by default

    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;

}; // class
///@}

///@} addtogroup block

} // namespace Kratos