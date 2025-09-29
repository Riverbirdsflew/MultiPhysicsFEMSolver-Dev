//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  Main authors:    whf
//                   



#if !defined(KRATOS_ISO)
#define  KRATOS_ISO

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "kinetics_model.h"
#include "../phase_transformation_type.h"


namespace Kratos
{

/**
 * Base class of trip.
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) Iso : public KineticsModel
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
     * Counted pointer of phase transformation model
     */
    KRATOS_CLASS_POINTER_DEFINITION(Iso);


    /**
     * Constructor.
     */
    Iso();

    /**
     * Destructor.
     */
    virtual ~Iso();

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this phase transformation model
     * @note implementation scheme:
     *      BasePhaseTransformationModel::Pointer p_clone(new BasePhaseTransformationModel());
     *      return p_clone;
     */
    virtual KineticsModel::Pointer Clone() const;


    /**
     * @return The working space dimension of the current phase transformation model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension();

     /**
      * Is called to check whether the provided material parameters in the Properties
      * match the requirements of current phase transformation model.
      * @param rMaterialProperties the current Properties to be validated against.
      * @return true, if parameters are correct; false, if parameters are insufficient / faulty
      * NOTE: this has to be implemented by each phase transformation model model. Returns false in base class since
      * no valid implementation is contained here.
      */
    virtual bool ValidateInput(const Properties& rMaterialProperties);


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
        buffer << "Iso" << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Iso";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
      rOStream << "Iso data";
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

}; /* Class PhaseTransformation */

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  Iso& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const Iso& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block


} /* namespace Kratos.*/
#endif /* KRATOS_AVRAMI  defined */
