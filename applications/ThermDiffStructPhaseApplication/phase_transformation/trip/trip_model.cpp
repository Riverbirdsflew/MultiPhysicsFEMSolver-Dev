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

#include "trip_model.h"

namespace Kratos
{

    /**
     * Constructor.
     */
    TripModel::TripModel()
    {
        this->mstat = Status::UnSet;
    }

    TripModel::~TripModel()
    {
        this->mstat = Status::UnSet;
    }
    
    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this trip model
     * @note implementation scheme:
     *      TripModel::Pointer p_clone(new TripModel());
     *      return p_clone;
     */
    TripModel::Pointer TripModel::Clone() const
    {
        KRATOS_ERROR <<  "Called the virtual function for clone from TripModel"<< std::endl;;
    }

    /**
     * @brief It creates a new trip model pointer
     * @param NewParameters The configuration parameters of the new trip model
     * @return a Pointer to the new trip model
     */
    
     TripModel::Pointer TripModel::Create(Kratos::Parameters NewParameters) const
     {
        const std::string& name = NewParameters["trip"].GetString();
        return KratosComponents<TripModel>::Get(name).Clone();
     }

     /**
     * @brief It creates a new trip model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new trip model
     * @param rProperties The properties of the material
     * @return a Pointer to the new trip model
     */
    TripModel::Pointer TripModel::Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const
    {
        return this->Create(NewParameters);
    }

    /**
     * @return The working space dimension of the current trip model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    TripModel::SizeType TripModel::WorkingSpaceDimension()
    {
        KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension from TripModel"<< std::endl;;
    }

    /**
     * @return The object status of the current kinetic model
     */
    const TripModel::Status TripModel::Stat()
    {
        /* return the flag of data status */
        return this->mstat;
    }

    int TripModel::Initialize(const Properties& rMaterialProperties,
                                    PhaseTransformationType rType)
    {
        KRATOS_ERROR << "Called the virtual function for initialize from TripModel" << std::endl;
        return 0;
    }

    int TripModel::Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR <<  "Called the virtual function for Check from TripModel"<< std::endl;;
        return 0;
    }

}
