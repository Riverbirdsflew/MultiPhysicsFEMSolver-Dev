//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  Main authors:    whf
//                   

#pragma once

#if !defined(KRATOS_PHASE_TRANSFORMATION_LAW)
#define  KRATOS_PHASE_TRANSFORMATION_LAW

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
#include "../phase_transformation/kinetics/kinetics_model.h"
#include "../phase_transformation/trip/trip_model.h"
#include "../phase_transformation/partition/partition_model.h"
#include "phase_transformation_type.h"


namespace Kratos
{

/**
 * Base class of phase transformation.
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) PhaseTransformationLaw : public Flags
{
public:


    /**
     * Type definitions
     * NOTE: geometries are assumed to be of type Node for all problems
     */
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node > GeometryType;
    typedef TripModel TTripType;
    typedef PartitionModel TPartitionType;
    typedef KineticsModel TKineticsType;

    /* define the status id */
	enum class Status : int {
		UnSet = 0,
		Set = 1,
	};


    /**
     * Counted pointer of phase transformation model
     */
    KRATOS_CLASS_POINTER_DEFINITION(PhaseTransformationLaw);

    

    struct ModelFeatures
    {
        KRATOS_CLASS_POINTER_DEFINITION(ModelFeatures);

        /**
         * Structure "Features" to be used by the element to get the phase transformation model characteristics*
         * its variables will be used to check phase transformation model and element compatibility
         */

        Flags           mOptions; // flags with the current phase transformation model characteristics
        SizeType        mStrainSize;
        SizeType        mSpaceDimension;

        /**
         * Constructor.
         */
        ModelFeatures() {};

        /**
         * Destructor.
         */
        ~ModelFeatures() {};

        // Set variables
        void SetOptions        (const Flags&  rOptions)        {mOptions=rOptions;};
        void SetStrainSize     (const SizeType StrainSize)     {mStrainSize=StrainSize;};
        void SetSpaceDimension (const SizeType SpaceDimension) {mSpaceDimension=SpaceDimension;};

        // Get variables
        const Flags& GetOptions () {return mOptions;};

        const SizeType& GetStrainSize()     {return mStrainSize;};
        const SizeType& GetSpaceDimension() {return mSpaceDimension;};
    };

    struct Parameters
    {
        KRATOS_CLASS_POINTER_DEFINITION(Parameters);

    /**
     * Structure "Parameters" to be used by the element to pass the parameters into the phase transformation model *

     * @param mOptions flags for the current phase transformation model Parameters (input data)

     * KINEMATIC PARAMETERS:

     *** NOTE: Pointers are used only to point to a certain variable, no "new" or "malloc" can be used for this Parameters ***

     * @param mpElementGeometry pointer to the element's geometry (input data)

     * MATERIAL PROPERTIES:
     * @param mpMaterialProperties pointer to the material's Properties object (input data)

     * PROCESS PROCESSINFO:
     * @param mpCurrentProcessInfo pointer to current ProcessInfo instance (input data)

     */

    private:

      /*** NOTE: Member Pointers are used only to point to a certain variable, no "new" or "malloc" can be used for this Parameters ***/

      Flags                                 mOptions;

      const Vector*                        mpShapeFunctionsValues;
      const ProcessInfo*                   mpCurrentProcessInfo;
      const Properties*                    mpMaterialProperties;
      const GeometryType*                  mpElementGeometry;

    public:


      /**
       * Constructor.
       */
      Parameters ()
      {
          //Initialize pointers to NULL
          mpShapeFunctionsValues=NULL;
          mpCurrentProcessInfo=NULL;
          mpMaterialProperties=NULL;
          mpElementGeometry=NULL;
      };


      /**
       * Constructor with Properties, Geometry and ProcessInfo
       */
      Parameters (
          const GeometryType& rElementGeometry,
          const Properties& rMaterialProperties,
          const ProcessInfo& rCurrentProcessInfo)
      :mpCurrentProcessInfo(&rCurrentProcessInfo)
      ,mpMaterialProperties(&rMaterialProperties)
      ,mpElementGeometry(&rElementGeometry)
      {
        mpShapeFunctionsValues = NULL;
      };

      /**
       * Copy Constructor.
       */
      Parameters (const Parameters & rNewParameters)
        :mOptions(rNewParameters.mOptions)
        ,mpShapeFunctionsValues(rNewParameters.mpShapeFunctionsValues)
        ,mpCurrentProcessInfo(rNewParameters.mpCurrentProcessInfo)
        ,mpMaterialProperties(rNewParameters.mpMaterialProperties)
        ,mpElementGeometry(rNewParameters.mpElementGeometry)
      {
      };

      /**
       * Destructor.
       */
      ~Parameters()
      {
      }

      /**
       * Verify Parameters
       */

      /**
       *Check currentprocessinfo, material properties and geometry
       */

      bool CheckInfoMaterialGeometry ()
      {
        if(!mpCurrentProcessInfo)
            KRATOS_ERROR << "CurrentProcessInfo NOT SET" << std::endl;

        if(!mpMaterialProperties)
            KRATOS_ERROR << "MaterialProperties NOT SET" << std::endl;

        if(!mpElementGeometry)
            KRATOS_ERROR << "ElementGeometry NOT SET" << std::endl;

        return 1;
      }



      /**
       * Public Methods to access variables of the struct class
       */

      /**
       * sets the variable or the pointer of a specified variable: assigns the direction of the pointer for the mpvariables, only non const values can be modified
       */

      void Set                             (Flags ThisFlag)                           {mOptions.Set(ThisFlag);};
      void Reset                           (Flags ThisFlag)                           {mOptions.Reset(ThisFlag);};

      void SetOptions                      (const Flags&  rOptions)                   {mOptions=rOptions;};
      void SetShapeFunctionsValues         (const Vector& rShapeFunctionsValues)      {mpShapeFunctionsValues=&rShapeFunctionsValues;};

      void SetProcessInfo                  (const ProcessInfo& rProcessInfo)          {mpCurrentProcessInfo =&rProcessInfo;};
      void SetMaterialProperties           (const Properties&  rMaterialProperties)   {mpMaterialProperties =&rMaterialProperties;};
      void SetElementGeometry              (const GeometryType& rElementGeometry)     {mpElementGeometry =&rElementGeometry;};

      /**
       * Returns the reference or the value of a specified variable: returns the value of the parameter, only non const values can be modified
       */
      Flags& GetOptions () {return mOptions;};

      const ProcessInfo& GetProcessInfo()
      {
          KRATOS_DEBUG_ERROR_IF_NOT(IsSetProcessInfo()) << "ProcessInfo is not set!" << std::endl;
          return *mpCurrentProcessInfo;
      }
      const Vector& GetShapeFunctionsValues()
      {
          KRATOS_DEBUG_ERROR_IF_NOT(IsSetShapeFunctionsValues()) << "ShapeFunctionsValues is not set!" << std::endl;
          return *mpShapeFunctionsValues;
      }
      const Properties& GetMaterialProperties()
      {
          KRATOS_DEBUG_ERROR_IF_NOT(IsSetMaterialProperties()) << "MaterialProperties is not set!" << std::endl;
          return *mpMaterialProperties;
      }
      const GeometryType& GetElementGeometry()
      {
          KRATOS_DEBUG_ERROR_IF_NOT(IsSetElementGeometry()) << "ElementGeometry is not set!" << std::endl;
          return *mpElementGeometry;
      }

      /**
       * Returns if the different components has been set
       */
      bool IsSetShapeFunctionsValues      () {return (mpShapeFunctionsValues != NULL);};
      bool IsSetProcessInfo               () {return (mpCurrentProcessInfo != NULL);};
      bool IsSetMaterialProperties        () {return (mpMaterialProperties != NULL);};
      bool IsSetElementGeometry           () {return (mpElementGeometry != NULL);};

    };// struct Parameters end

    /**
     * Constructor.
     */
    PhaseTransformationLaw();

    /**
     * Constructor with model information.
     */
    PhaseTransformationLaw(
        PhaseTransformationType type,
        TKineticsType::Pointer pKineticsModel,
        TPartitionType::Pointer pPartitionModel,
        TTripType::Pointer pTripModel);

    /**
     * Copy Constructor.
     */
    PhaseTransformationLaw(const PhaseTransformationLaw& other);
        

    /**
     * Destructor.
     */
    ~PhaseTransformationLaw() override = default;

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this phase transformation model
     * @note implementation scheme:
     *      PhaseTransformationLaw::Pointer p_clone(new PhaseTransformationLaw());
     *      return p_clone;
     */
    virtual PhaseTransformationLaw::Pointer Clone() const;

    /**
     * @brief It creates a new phase transformation model pointer
     * @param NewParameters The configuration parameters of the new phase transformation model
     * @return a Pointer to the new phase transformation model
     */
    virtual Pointer Create(Kratos::Parameters NewParameters) const;

    /**
     * @brief It creates a new phase transformation model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new phase transformation model
     * @param rProperties The properties of the material
     * @return a Pointer to the new phase transformation model
     */
    virtual Pointer Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const;

    /**
     * @return The working space dimension of the current phase transformation model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension();


    /**
     * Get internal members
     */
    virtual const Status Stat();
    virtual const Variable<double>* GetLatentPtr();
	virtual const Variable<double>* GetTransStrainPtr();
	virtual const Variable<double>* GetUpperTPtr();
	virtual const Variable<double>* GetLowerTPtr();
	virtual const Variable<double>* GetMaxMassFracPtr();
    virtual const TKineticsType::Pointer GetKineticsModel();
    virtual const TPartitionType::Pointer GetPartitionModel();
    virtual const TTripType::Pointer GetTripModel();

    virtual PhaseTransformationType GetPhaseTranType() const
    {
        return this->mPhaseTranType;
    }
    // 返回相变类型的字符串
    std::string GetPhaseTranTypeName() const
    {
        switch (mPhaseTranType) 
        {
            case PhaseTransformationType::AUSTENIZATION:   return "AUSTENIZATION";
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:   return "EUTECTOID_DECOMPOSITION";
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:     return "BAINITIC_TRANSFORMATION";
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:  return "MARTENSITIC_TRANSFORMATION";
            default:                         return "Unknown Transformation";
        }
    }

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the phase transformation model
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual int Initialize(const Properties& rMaterialProperties, PhaseTransformationType rType);

    /**
     * @brief Returns the mean value integral value of unit latent of phase transformation
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rCurrentProcessInfo the current ProcessInfo instance
     * @return the mean value of the latent heat of phase transformation
     */
    virtual double GetPhaseTransLatent(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

    virtual double GetPhaseTransStrain(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

    virtual double GetPhaseTransUpperT(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

    virtual double GetPhaseTransLowerT(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

    virtual double GetMaxTransMassFrc(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      SizeType NodeId);

    virtual int GetPartitionInfo(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double rMassFrac, double *rMassTransFracInc);

    virtual int GetPhaseTransMassFrcInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                        const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                        SizeType NodeId, double rIncuFrac, double rMassFrac, double rSatuFrac,
                                        double &rIncuFracInc, double &rMassFracInc);

    virtual double GetEqTripInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                SizeType NodeId, double rMassFrac, double rMassFracInc);


    /**
     * @brief Initialize the material response,  called by the element in InitializeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */
    void InitializeMaterialResponse (Parameters& rValues);

    /**
     * @brief Finalize the material response,  called by the element in FinalizeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */
    void FinalizeMaterialResponse (Parameters& rValues);

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rModelFeatures
     */
    virtual void GetModelFeatures(ModelFeatures& rModelFeatures);

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
    std::string Info() const override
    {
        return "PhaseTransformationLaw: " + this->GetPhaseTranTypeName();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "PhaseTransformationLaw has no data" << std::endl;
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
        PhaseTransformationType mPhaseTranType;
        Variable<double>* mLatent;
        Variable<double>* mTransStrain;
        Variable<double>* mUpperT;
        Variable<double>* mLowerT;
        Variable<double>* mMaxMassFrc;
        TPartitionType::Pointer mpPartition;
        TKineticsType::Pointer mpKinetics;
        TTripType::Pointer mpTrip;

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

}; /* Class Phasetransformation */

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  PhaseTransformationLaw& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const PhaseTransformationLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

template class KratosComponents<PhaseTransformationLaw >;
template class KratosComponents< Variable<PhaseTransformationLaw::Pointer> >;

void AddKratosComponent(std::string const& Name, PhaseTransformationLaw const& ThisComponent);
void AddKratosComponent(std::string const& Name, Variable<PhaseTransformationLaw::Pointer> const& ThisComponent);

/**
 * Definition of phase transformation model variable
 */

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(PhaseTransformationLaw::Pointer, PHASE_TRANSFORMATION_LAW)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT


} /* namespace Kratos.*/
#endif /* KRATOS_PHASE_TRANSFORMATION_LAW  defined */
