//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  Main authors:    whf
//                   



#if !defined(KRATOS_KME)
#define  KRATOS_KME

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "kinetics_model.h"
#include "../phase_transformation_law.h"

namespace Kratos
{

/**
 * Base class of trip.
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) Kme : public KineticsModel
{
public:

    /**
     * Type definitions
     * NOTE: geometries are assumed to be of type Node for all problems
     */
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node > GeometryType;


    /**
     * Counted pointer of phase transition model
     */
    KRATOS_CLASS_POINTER_DEFINITION(Kme);


    /**
     * Constructor.
     */
    Kme();

    /**
     * Destructor.
     */
    virtual ~Kme();

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this phase transition model
     * @note implementation scheme:
     *      BasePhaseTransitionModel::Pointer p_clone(new BasePhaseTransitionModel());
     *      return p_clone;
     */
    virtual KineticsModel::Pointer Clone() const;


    /**
     * @return The working space dimension of the current phase transition model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension();

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the trip model
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual int Initialize(const Properties& rMaterialProperties, PhaseTransformationType rType);

    virtual int GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
                                        const GeometryType &rElementGeometry,
                                        const Vector &rShapeFunctionsValues,
                                        const ProcessInfo &rCurrentProcessInfo,
                                        SizeType NodeId,
                                        double rIncuFrac,
                                        double rMassFrac,
                                        double rSatuFrac,
                                        double &rIncuFracInc,
                                        double &rMassFracInc);

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
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Kme" << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Kme";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
      rOStream << "Kme data";
    }


    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    Variable<double>* mPara[2]; 
    PhaseTransformationType mType;

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
                                  Kme& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const Kme& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block


} /* namespace Kratos.*/
#endif /* KRATOS_KME  defined */
