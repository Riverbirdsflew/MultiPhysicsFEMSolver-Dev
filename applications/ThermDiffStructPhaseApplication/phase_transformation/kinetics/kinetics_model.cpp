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

#include "kinetics_model.h"
#include "phase_transformation_type.h"

namespace Kratos
{

    /**
     * Constructor.
     */
    KineticsModel::KineticsModel()
    {
        this->mstat = Status::UnSet;
    }

    KineticsModel::~KineticsModel()
    {
        this->mstat = Status::UnSet;
    }
    
    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this kinetics model
     * @note implementation scheme:
     *      KineticsModel::Pointer p_clone(new KineticsModel());
     *      return p_clone;
     */
    KineticsModel::Pointer KineticsModel::Clone() const
    {
        KRATOS_ERROR <<  "Called the virtual function for clone from KineticsModel"<< std::endl;;
    }

    /**
     * @brief It creates a new kinetics model pointer
     * @param NewParameters The configuration parameters of the new kinetics model
     * @return a Pointer to the new kinetics model
     */
    
     KineticsModel::Pointer KineticsModel::Create(Kratos::Parameters NewParameters) const
     {
        const std::string& name = NewParameters["kinetics"].GetString();
        return KratosComponents<KineticsModel>::Get(name).Clone();
     }

     /**
     * @brief It creates a new kinetics model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new kinetics model
     * @param rProperties The properties of the material
     * @return a Pointer to the new kinetics model
     */
    KineticsModel::Pointer KineticsModel::Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const
    {
        return this->Create(NewParameters);
    }

    /**
     * @return The working space dimension of the current kinetics model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    KineticsModel::SizeType KineticsModel::WorkingSpaceDimension()
    {
        KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension from KineticsModel"<< std::endl;;
    }

    /**
     * @return The object status of the current kinetic model
     */
    const KineticsModel::Status KineticsModel::Stat()
    {
        /* return the flag of data status */
        return this->mstat;
    }


    int KineticsModel::Initialize(const Properties& rMaterialProperties,
                                    PhaseTransformationType rType)
    {
        KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension from KineticsModel"<< std::endl;;
        return 0;
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int KineticsModel::Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR <<  "Called the virtual function for Check from KineticsModel"<< std::endl;;
        return 0; // Placeholder return value, should be replaced with actual implementation
    }

}
