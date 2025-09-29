//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  Main authors:    whf
//                   



#if !defined(KRATOS_PARTITION_MODEL)
#define  KRATOS_PARTITION_MODEL

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "containers/data_value_container.h"
#include "containers/flags.h"
#include "includes/initial_state.h"
#include "../phase_transformation_type.h"


namespace Kratos
{

/**
 * Base class of trip.
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) PartitionModel
{
public:

    /**
     * Type definitions
     * NOTE: geometries are assumed to be of type Node for all problems
     */
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node > GeometryType;

    /* define the status id */
	enum class Status : int {
		UnSet = 0,
		Set = 1,
	};

    /**
     * Counted pointer of phase transition model
     */
    KRATOS_CLASS_POINTER_DEFINITION(PartitionModel);


    /**
     * Constructor.
     */
    PartitionModel();

    /**
     * Destructor.
     */
    virtual ~PartitionModel();

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this phase transition model
     */
    virtual PartitionModel::Pointer Clone() const;

    /**
     * @brief It creates a new phase transition model pointer
     * @param NewParameters The configuration parameters of the new phase transition model
     * @return a Pointer to the new phase transition model
     */
    virtual Pointer Create(Kratos::Parameters NewParameters) const;

    /**
     * @brief It creates a new phase transition model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new phase transition model
     * @param rProperties The properties of the material
     * @return a Pointer to the new phase transition model
     */
    virtual Pointer Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const;

    /**
     * @return The working space dimension of the current phase transition model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension();

    /* return the internal variables */
	virtual const Status Stat();

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the trip model
     */
    virtual int Initialize(const Properties& rMaterialProperties, PhaseTransformationType rType);

    virtual int GetPartitionInfo(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                 const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                 SizeType NodeId, double rMassFrac, double *rMassTransFracInc);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    virtual int Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const;


    ///@}
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PartitionModel";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PartitionModel" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "PartitionModel has no data" << std::endl;
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

        Status mstat;

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


    ///@}
    ///@name Member Variables
    ///@{

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

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class PhaseTransition */

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  PartitionModel& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const PartitionModel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

template class KratosComponents<PartitionModel >;
template class KratosComponents< Variable<PartitionModel::Pointer> >;

void AddKratosComponent(std::string const& Name, PartitionModel const& ThisComponent);
void AddKratosComponent(std::string const& Name, Variable<PartitionModel::Pointer> const& ThisComponent);

/**
 * Definition of trip model variable
 */

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(PartitionModel::Pointer, PARTITION_MODEL)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT


} /* namespace Kratos.*/
#endif /* KRATOS_BASE_PHASE_TRANSITION_MODEL  defined */
