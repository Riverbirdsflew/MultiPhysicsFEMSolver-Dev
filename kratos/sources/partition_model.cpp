//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    whf
//

#include "includes/partition_model.h"

namespace Kratos
{

    /**
     * Constructor.
     */
    PartitionModel::PartitionModel() : Flags()
    {
    }
    
    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this partition model
     * @note implementation scheme:
     *      PartitionModel::Pointer p_clone(new PartitionModel());
     *      return p_clone;
     */
    PartitionModel::Pointer PartitionModel::Clone() const
    {
        KRATOS_ERROR <<  "Called the virtual function for clone from PartitionModel"<< std::endl;;
    }

    /**
     * @brief It creates a new partition model pointer
     * @param NewParameters The configuration parameters of the new partition model
     * @return a Pointer to the new partition model
     */
    
     PartitionModel::Pointer PartitionModel::Create(Kratos::Parameters NewParameters) const
     {
        const std::string& name = NewParameters["name"].GetString();
        return KratosComponents<PartitionModel>::Get(name).Clone();
     }

     /**
     * @brief It creates a new partition model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new partition model
     * @param rProperties The properties of the material
     * @return a Pointer to the new partition model
     */
    PartitionModel::Pointer PartitionModel::Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const
    {
        return this->Create(NewParameters);
    }

    /**
     * @return The working space dimension of the current partition model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    PartitionModel::SizeType PartitionModel::WorkingSpaceDimension()
    {
        KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension from PartitionModel"<< std::endl;
    }


    int PartitionModel::Initialize(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const Vector& rShapeFunctionsValues)
    {
        KRATOS_ERROR <<  "Called the virtual function for Initialize from PartitionModel"<< std::endl;
        return -1; // Placeholder return value, should be replaced with actual implementation
    }

    int PartitionModel::GetPartitionInfo(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR <<  "Called the virtual function for GetPartitionInfo from PartitionModel"<< std::endl;
        return -1; // Placeholder return value, should be replaced with actual implementation
    }

    int PartitionModel::Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR <<  "Called the virtual function for Check from PartitionModel"<< std::endl;
        return -1; // Return an error code
    }

}
