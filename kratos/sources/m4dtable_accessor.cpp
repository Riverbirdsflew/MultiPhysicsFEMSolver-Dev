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

// System includes

// External includes

// Project includes
#include "includes/m4dtable_accessor.h"
#include "includes/properties.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetValue(
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetValueFromM4dTable(mDim, mpInputVariable, rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetValue(
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const Node& rNode,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetValueFromM4dTable(mDim, mpInputVariable, rVariable, rProperties, rNode, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetDerivative(
    const Variable<double>& rDerVariable,
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetDerivativeFromM4dTable(mDim, mpInputVariable, rDerVariable, rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetDerivative(
    const Variable<double>& rDerVariable,
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const Node& rNode,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetDerivativeFromM4dTable(mDim, mpInputVariable, rDerVariable, rVariable, rProperties, rNode, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetIntegration(
    const Variable<double>& rIntVariable,
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetIntegrationFromM4dTable(mDim, mpInputVariable, rIntVariable, rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetIntegration(
    const Variable<double>& rIntVariable,
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const Node& rNode,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetIntegrationFromM4dTable(mDim, mpInputVariable, rIntVariable, rVariable, rProperties, rNode, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetValueFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const
{
    M4dTable::ValueListType independent_at_gauss;
    SizeType m, i;

    // init the independent_at_gauss
    for (m = 0; m < 4; m++) {
        independent_at_gauss[m] = 0.0;
    }

    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            for (i = 0; i < rShapeFunctionVector.size(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].SolutionStepsDataHas(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
                const double nodal_value = rGeometry[i].FastGetSolutionStepValue(rp);
                independent_at_gauss[m] += nodal_value * rShapeFunctionVector[i];
            }
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            for (i = 0; i < rShapeFunctionVector.size(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].Has(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
                const double nodal_value = rGeometry[i].GetValue(rp);
                independent_at_gauss[m] += nodal_value * rShapeFunctionVector[i];
            }
        }
    } else if (mInputVariableType == Globals::DataLocation::Element) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry.Has(rp)) << "The Variable " << rp.Name() << " is not available at the Geometry to retrieve M4dTable values." << std::endl;
            independent_at_gauss[m] = rGeometry.GetValue(rp);
        }
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical, nodal_non_historical and elemental_non_historical" << std::endl;
    }

    // Retrieve the dependent variable from the table
    if (rDim == 1) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(independent_at_gauss);
    } else if (rDim == 2) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(independent_at_gauss);
    } else if (rDim == 3) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(independent_at_gauss);
    } else if (rDim == 4) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& rIndenpendentVariable4 = *rpIndependentVariable[3];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable4, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(independent_at_gauss);
    } else {
        return 0.0;
    }    
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetValueFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo) const
{
    M4dTable::ValueListType nodal_value;
    SizeType m;

    // init the nodal_value
    for (m = 0; m < 4; m++) {
        nodal_value[m] = 0.0;
    }

    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];
            KRATOS_DEBUG_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
            nodal_value[m] = rNode.FastGetSolutionStepValue(rp);
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];
            KRATOS_DEBUG_ERROR_IF_NOT(rNode.Has(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
            nodal_value[m] = rNode.GetValue(rp);
        }
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical and nodal_non_historical" << std::endl;
    }

    // Retrieve the dependent variable from the table
    if (rDim == 1) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(nodal_value);
    } else if (rDim == 2) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(nodal_value);
    } else if (rDim == 3) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(nodal_value);
    } else if (rDim == 4) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& rIndenpendentVariable4 = *rpIndependentVariable[3];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable4, rDependentVariable);
        return r_m4dtable.GetInterpolateValue(nodal_value);
    } else {
        return 0.0;
    }    
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetDerivativeFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDerIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const
{
    M4dTable::ValueListType independent_at_gauss;
    SizeType m, n, i;

    // init the independent_at_gauss and check the derivative independent variable
    for (m = 0, n = 0; m < 4; m++) {
        independent_at_gauss[m] = 0.0;

        if (rpIndependentVariable[m] == &rDerIndependentVariable)
            n = m + 1;
    }

    // check the dimension of derivative independent variable
    if (n == 0)  
        return 0.0;
    else
        n--;

    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            for (i = 0; i < rShapeFunctionVector.size(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].SolutionStepsDataHas(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
                const double nodal_value = rGeometry[i].FastGetSolutionStepValue(rp);
                independent_at_gauss[m] += nodal_value * rShapeFunctionVector[i];
            }
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            for (i = 0; i < rShapeFunctionVector.size(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].Has(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
                const double nodal_value = rGeometry[i].GetValue(rp);
                independent_at_gauss[m] += nodal_value * rShapeFunctionVector[i];
            }
        }
    } else if (mInputVariableType == Globals::DataLocation::Element) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry.Has(rp)) << "The Variable " << rp.Name() << " is not available at the Geometry to retrieve M4dTable values." << std::endl;
            independent_at_gauss[m] = rGeometry.GetValue(rp);
        }
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical, nodal_non_historical and elemental_non_historical" << std::endl;
    }

    // Retrieve the dependent derivative variable from the table
    if (rDim == 1) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, independent_at_gauss);
    } else if (rDim == 2) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, independent_at_gauss);
    } else if (rDim == 3) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, independent_at_gauss);
    } else if (rDim == 4) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& rIndenpendentVariable4 = *rpIndependentVariable[3];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable4, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, independent_at_gauss);
    } else {
        return 0.0;
    } 
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetDerivativeFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rDerIndependentVariable,
        const Variable<double>& rDependentVariable,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo) const
{
    M4dTable::ValueListType nodal_value;
    SizeType m, n;

    // init the nodal_value and check the derivative independent variable
    for (m = 0, n = 0; m < 4; m++) {
        nodal_value[m] = 0.0;

        if (rpIndependentVariable[m] == &rDerIndependentVariable)
            n = m + 1;
    }

    // check the dimension of derivative independent variable
    if (n == 0)  
        return 0.0;
    else
        n--;

    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];
            KRATOS_DEBUG_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
            nodal_value[m] = rNode.FastGetSolutionStepValue(rp);
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];
            KRATOS_DEBUG_ERROR_IF_NOT(rNode.Has(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
            nodal_value[m] = rNode.GetValue(rp);
        }
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical and nodal_non_historical" << std::endl;
    }

    // Retrieve the dependent derivative variable from the table
    if (rDim == 1) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, nodal_value);
    } else if (rDim == 2) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, nodal_value);
    } else if (rDim == 3) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, nodal_value);
    } else if (rDim == 4) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& rIndenpendentVariable4 = *rpIndependentVariable[3];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable4, rDependentVariable);
        return r_m4dtable.GetDerivativeValue(n, nodal_value);
    } else {
        return 0.0;
    } 
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetIntegrationFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rIntIndependentVariable,
        const Variable<double> &rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const
{
    M4dTable::ValueListType independent_at_gauss;
    M4dTable::ValueListType init_independent_at_gauss;
    SizeType m, n, i;

    // init the independent_at_gauss and check the derivative independent variable
    for (m = 0, n = 0; m < 4; m++) {
        independent_at_gauss[m] = 0.0;
        init_independent_at_gauss[m] = 0.0;

        if (rpIndependentVariable[m] == &rIntIndependentVariable)
            n = m + 1;
    }

    // check the dimension of derivative independent variable
    if (n == 0)  
        return 0.0;
    else
        n--;

    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            for (i = 0; i < rShapeFunctionVector.size(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].SolutionStepsDataHas(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
                const double nodal_value = rGeometry[i].FastGetSolutionStepValue(rp);
                independent_at_gauss[m] += nodal_value * rShapeFunctionVector[i];
            }
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            for (i = 0; i < rShapeFunctionVector.size(); ++i) {
                KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].Has(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
                const double nodal_value = rGeometry[i].GetValue(rp);
                independent_at_gauss[m] += nodal_value * rShapeFunctionVector[i];
            }
        }
    } else if (mInputVariableType == Globals::DataLocation::Element) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];

            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry.Has(rp)) << "The Variable " << rp.Name() << " is not available at the Geometry to retrieve M4dTable values." << std::endl;
            independent_at_gauss[m] = rGeometry.GetValue(rp);
        }
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical, nodal_non_historical and elemental_non_historical" << std::endl;
    }

    for(m = 0; m < rDim; m++)
    {
        if(m != n)  
            init_independent_at_gauss[m] =  independent_at_gauss[m];
    }

    // Retrieve the dependent derivative variable from the table
    if (rDim == 1) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_independent_at_gauss, independent_at_gauss);
    } else if (rDim == 2) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_independent_at_gauss, independent_at_gauss);
    } else if (rDim == 3) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_independent_at_gauss, independent_at_gauss);
    } else if (rDim == 4) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& rIndenpendentVariable4 = *rpIndependentVariable[3];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable4, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_independent_at_gauss, independent_at_gauss);
    } else {
        return 0.0;
    } 
}

/***********************************************************************************/
/***********************************************************************************/

double M4dTableAccessor::GetIntegrationFromM4dTable(
        const SizeType& rDim,
        const VariableListType& rpIndependentVariable,
        const Variable<double>& rIntIndependentVariable,
        const Variable<double> &rDependentVariable,
        const Properties& rProperties,
        const Node& rNode,
        const ProcessInfo& rProcessInfo) const
{
    M4dTable::ValueListType nodal_value;
    M4dTable::ValueListType init_nodal_value;
    SizeType m, n;

    // init the nodal_value and check the derivative independent variable
    for (m = 0, n = 0; m < 4; m++) {
        nodal_value[m] = 0.0;
        init_nodal_value[m] = 0.0;

        if (rpIndependentVariable[m] == &rIntIndependentVariable)
            n = m + 1;
    }

    // check the dimension of derivative independent variable
    if (n == 0)  
        return 0.0;
    else
        n--;

    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];
            KRATOS_DEBUG_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
            nodal_value[m] = rNode.FastGetSolutionStepValue(rp);
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for(m = 0; m < rDim; m++) {
            auto& rp = *rpIndependentVariable[m];
            KRATOS_DEBUG_ERROR_IF_NOT(rNode.Has(rp)) << "The Variable " << rp.Name() << " is not available at the nodes of the Geometry to retrieve M4dTable values." << std::endl;
            nodal_value[m] = rNode.GetValue(rp);
        }
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical and nodal_non_historical" << std::endl;
    }

    for(m = 0; m < rDim; m++)
    {
        if(m != n)  
            init_nodal_value[m] =  nodal_value[m];
    }

    // Retrieve the dependent derivative variable from the table
    if (rDim == 1) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_nodal_value, nodal_value);
    } else if (rDim == 2) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable1, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_nodal_value, nodal_value);
    } else if (rDim == 3) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable1, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_nodal_value, nodal_value);
    } else if (rDim == 4) {
        const auto& rIndenpendentVariable1 = *rpIndependentVariable[0];
        const auto& rIndenpendentVariable2 = *rpIndependentVariable[1];
        const auto& rIndenpendentVariable3 = *rpIndependentVariable[2];
        const auto& rIndenpendentVariable4 = *rpIndependentVariable[3];
        const auto& r_m4dtable = rProperties.GetM4dTable(rDim, rIndenpendentVariable1, rIndenpendentVariable2, rIndenpendentVariable3, rIndenpendentVariable4, rDependentVariable);
        return r_m4dtable.GetIntegrationValue(init_nodal_value, nodal_value);
    } else {
        return 0.0;
    } 
}

/***********************************************************************************/
/***********************************************************************************/

void M4dTableAccessor::save(Serializer& rSerializer) const
{
    rSerializer.save("InputVariableDimension", mDim);
  
    for (SizeType m = 0; m < mDim; m++)
        rSerializer.save("InputVariable", mpInputVariable[m]->Name());

    // // we must do the int cast to be able to compile
    rSerializer.save("InputVariableType", static_cast<int>(mInputVariableType)); 
}

void M4dTableAccessor::load(Serializer& rSerializer)
{
    std::string variable_name;
    int enum_value;

    rSerializer.load("InputVariableDimension", mDim);

    for (SizeType m = 0; m < mDim; m++)
    {
        rSerializer.load("InputVariable", variable_name);
        mpInputVariable[m] = static_cast<Variable<double> *>(KratosComponents<VariableData>::pGet(variable_name));
    }

    // // we must do the int cast to be able to compile   
    rSerializer.load("InputVariableType", enum_value);
    mInputVariableType = static_cast<Globals::DataLocation>(enum_value);
}

/***********************************************************************************/
/***********************************************************************************/

Accessor::UniquePointer M4dTableAccessor::Clone() const
{
    return Kratos::make_unique<M4dTableAccessor>(*this);
}

} // namespace Kratos
