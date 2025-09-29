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

#include "includes/creep_model.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    CreepModel::CreepModel() : Flags()
    {
    }
    
    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this creep model
     * @note implementation scheme:
     *      CreepModel::Pointer p_clone(new CreepModel());
     *      return p_clone;
     */
    CreepModel::Pointer CreepModel::Clone() const
    {
        KRATOS_ERROR <<  "Called the virtual function for clone from CreepModel"<< std::endl;;
    }

    /**
     * @brief It creates a new creep model pointer
     * @param NewParameters The configuration parameters of the new creep model
     * @return a Pointer to the new creep model
     */
    
     CreepModel::Pointer CreepModel::Create(Kratos::Parameters NewParameters) const
     {
        const std::string& name = NewParameters["name"].GetString();
        return KratosComponents<CreepModel>::Get(name).Clone();
     }

     /**
     * @brief It creates a new creep model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new creep model
     * @param rProperties The properties of the material
     * @return a Pointer to the new creep model
     */
    CreepModel::Pointer CreepModel::Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const
    {
        return this->Create(NewParameters);
    }

    /**
     * @return The working space dimension of the current creep model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    CreepModel::SizeType CreepModel::WorkingSpaceDimension()
    {
        KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension from CreepModel"<< std::endl;
    }


    int CreepModel::Initialize(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const Vector& rShapeFunctionsValues)
    {
        KRATOS_ERROR <<  "Called the virtual function for init from CreepModel"<< std::endl;
        return -1; // Return an error code
    }

    int CreepModel::GetCreepParameters(const Properties& rMaterialProperties,
                           const GeometryType& rElementGeometry,
                           const ProcessInfo& rCurrentProcessInfo,
                                unsigned int *pnm, double *pdata)
    {
        KRATOS_ERROR <<  "Called the virtual function for GetCreepParameters from CreepModel"<< std::endl;
        return -1; // Return an error code
    }

    double CreepModel::GetEqCreepRate(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR <<  "Called the virtual function for GetEqCreepRate from CreepModel"<< std::endl;
        return -1.0; // Return an error code
    }
    double CreepModel::GetEqCreepRateMean(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR <<  "Called the virtual function for GetEqCreepRateMean from CreepModel"<< std::endl;
        return -1.0; // Return an error code
    }

    int CreepModel::Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR <<  "Called the virtual function for Check from CreepModel"<< std::endl;
        return -1; // Return an error code
    }
}
