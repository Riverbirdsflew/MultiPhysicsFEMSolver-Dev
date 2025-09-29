//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//		Marx Xu

// System includes

// External includes

// Project includes
#include "utilities/read_and_set_accessors_utilities.h"

namespace Kratos {

/***********************************************************************************/
/***********************************************************************************/

void ReadAndSetAccessorsUtilities::ReadAndSetAccessors(
    const Parameters MaterialData,
    Properties& rProperty
    )
{
    if (MaterialData.Has("accessors")) {
        Parameters accessors = MaterialData["accessors"];

        // Loop over the accessors list
        for (auto iter = accessors.begin(); iter != accessors.end(); ++iter) {
            auto accessor_param = accessors.GetValue(iter.name());
            
            // Table Accessor
            if (accessor_param["accessor_type"].GetString() == "table_accessor") {
                // Independent Variable
                std::string input_var_name = accessor_param["properties"]["table_input_variable"].GetString();
                Variable<double> *p_input_var = static_cast<Variable<double> *>(KratosComponents<VariableData>::pGet(input_var_name));

                // Dependent Variable
                std::string output_var_name = accessor_param["properties"]["table_output_variable"].GetString();
                const auto& r_output_var  = KratosComponents<Variable<double>>().Get(output_var_name);

                // We set the variable type of the input variable (node_historical, node_non_historical and element)
                std::string input_var_type = accessor_param["properties"].Has("table_input_variable_type") ? accessor_param["properties"]["table_input_variable_type"].GetString() : "node_historical";

                KRATOS_ERROR_IF(rProperty.HasAccessor(r_output_var)) << "You are trying to add an TableAccessor between " << input_var_name << " and " << output_var_name << " which already exists..." << std::endl;

                rProperty.SetAccessor(r_output_var, (TableAccessor(*p_input_var, input_var_type)).Clone());
            } else if (accessor_param["accessor_type"].GetString() == "m4dtable_accessor") {
	// Dimension of Independent Variable
	M4dTableAccessor::SizeType rdim = accessor_param["properties"]["m4dtable_input_dimension"].GetInt();

	// Independent Variable
	Variable<double>* p_input_var[4];

	for (M4dTableAccessor::SizeType i = 0; i < rdim; i++)
	{
	    std::string input_var_tag = "m4dtable_input_variable_";
	    input_var_tag += std::to_string(i + 1);
	    std::string input_var_name = accessor_param["properties"][input_var_tag.c_str()].GetString();
                    p_input_var[i] = static_cast<Variable<double>*>(KratosComponents<VariableData>::pGet(input_var_name));
	}

                // Dependent Variable
                std::string output_var_name = accessor_param["properties"]["m4dtable_output_variable"].GetString();
                const auto& r_output_var  = KratosComponents<Variable<double>>().Get(output_var_name);

                // We set the variable type of the input variable (node_historical, node_non_historical and element)
                std::string input_var_type = accessor_param["properties"].Has("m4dtable_input_variable_type") ? accessor_param["properties"]["m4dtable_input_variable_type"].GetString() : "node_historical";

                KRATOS_ERROR_IF(rProperty.HasAccessor(r_output_var)) << "You are trying to add an M4dTableAccessor between " << p_input_var[0]->Name() << " and " << output_var_name << " which already exists..." << std::endl;

                rProperty.SetAccessor(r_output_var, (M4dTableAccessor(rdim, p_input_var, input_var_type)).Clone());
            } else {
                // No more accessors implemented currently
                KRATOS_ERROR << "This Accessor type is not available, only TableAccessor is ready for now" << std::endl;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos.
