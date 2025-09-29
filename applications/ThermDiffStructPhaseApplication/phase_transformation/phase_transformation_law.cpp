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

#include "phase_transformation_law.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    PhaseTransformationLaw::PhaseTransformationLaw()
    {
        this->mstat = Status::UnSet;
        this->mPhaseTranType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mLatent = nullptr;
        this->mTransStrain = nullptr;
        this->mUpperT = nullptr;
        this->mLowerT = nullptr;
        this->mMaxMassFrc = nullptr;
        this->mpKinetics = nullptr;
        this->mpPartition = nullptr;
        this->mpTrip = nullptr;
    }

    /**
     * Constructor with model information.
     */
    PhaseTransformationLaw::PhaseTransformationLaw(
        PhaseTransformationType type,
        TKineticsType::Pointer pKineticsModel,
        TPartitionType::Pointer pPartitionModel,
        TTripType::Pointer pTripModel)
    {
        this->mstat = Status::UnSet;
        this->mPhaseTranType = type;
        this->mLatent = nullptr;
        this->mTransStrain = nullptr;
        this->mUpperT = nullptr;
        this->mLowerT = nullptr;
        this->mMaxMassFrc = nullptr;
        this->mpKinetics = pKineticsModel;
        this->mpPartition = pPartitionModel;
        this->mpTrip = pTripModel;
    }

    /**
     * Copy Constructor.
     */
    PhaseTransformationLaw::PhaseTransformationLaw(const PhaseTransformationLaw& other)
    {
        this->mstat = other.mstat;
        this->mPhaseTranType = other.mPhaseTranType;
        this->mLatent = other.mLatent;
        this->mTransStrain = other.mTransStrain;
        this->mUpperT = other.mUpperT;
        this->mLowerT = other.mLowerT;
        this->mMaxMassFrc = other.mMaxMassFrc;
        this->mpKinetics = other.mpKinetics;
        this->mpPartition = other.mpPartition;
        this->mpTrip = other.mpTrip;
    }
    /**
     * Destructor.
     */
    PhaseTransformationLaw::~PhaseTransformationLaw()
    {
        this->mstat = Status::UnSet;
        this->mPhaseTranType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mLatent = nullptr;
        this->mTransStrain = nullptr;
        this->mUpperT = nullptr;
        this->mLowerT = nullptr;
        this->mMaxMassFrc = nullptr;
        this->mpKinetics = nullptr;
        this->mpPartition = nullptr;
        this->mpTrip = nullptr;
    }

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this phase transformation model
     * @note implementation scheme:
     *      PhaseTransformationLaw::Pointer p_clone(new PhaseTransformationLaw());
     *      return p_clone;
     */
    PhaseTransformationLaw::Pointer PhaseTransformationLaw::Clone() const
    {
        return Kratos::make_shared<PhaseTransformationLaw>(*this);
    }

    /**
     * @brief It creates a new phase transformation model pointer
     * @param NewParameters The configuration parameters of the new phase transformation model
     * @return a Pointer to the new phase transformation model
     */

    PhaseTransformationLaw::Pointer PhaseTransformationLaw::Create(Kratos::Parameters NewParameters) const
    {
        PhaseTransformationType type = PhaseTransformationType::UNDEFINED_TRANSFORMATION;

        if (NewParameters.Has("type"))
        {
            const std::string &type_name = NewParameters["type"].GetString();
            type = GetPhaseTransTypeFromString(type_name);
        }

        return Kratos::make_shared<PhaseTransformationLaw>(type);
    }

     /**
     * @brief It creates a new phase transformation model pointer (version with properties)
     * @param NewParameters The configuration parameters of the new phase transformation model
     * @param rProperties The properties of the material
     * @return a Pointer to the new phase transformation model
     */
    PhaseTransformationLaw::Pointer PhaseTransformationLaw::Create(
        Kratos::Parameters NewParameters,
        const Properties& rProperties
        ) const
    {
        return this->Create(NewParameters);
    }

    /**
     * @return The working space dimension of the current phase transformation model
     * @note This function HAS TO BE IMPLEMENTED by any derived class
     */
    PhaseTransformationLaw::SizeType PhaseTransformationLaw::WorkingSpaceDimension()
    {
        KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension from PhaseTransformationLaw"<< std::endl;;
    }

    /**
     * Get internal members
     */
    const PhaseTransformationLaw::Status PhaseTransformationLaw::Stat()
    {
        /* return the flag of data status */
        return this->mstat;
    }

    const Variable<double>* PhaseTransformationLaw::GetLatentPtr()
    {
        return (this->mstat == Status::Set ? this->mLatent : nullptr);
    }
    const Variable<double>* PhaseTransformationLaw::GetTransStrainPtr()
    {
        return (this->mstat == Status::Set ? this->mTransStrain : nullptr);
    }
    const Variable<double>* PhaseTransformationLaw::GetUpperTPtr()
    {
        return (this->mstat == Status::Set ? this->mUpperT : nullptr);
    }
    const Variable<double>* PhaseTransformationLaw::GetLowerTPtr()
    {
        return (this->mstat == Status::Set ? this->mLowerT : nullptr);
    }
    const Variable<double>* PhaseTransformationLaw::GetMaxMassFracPtr()
    {
        return (this->mstat == Status::Set ? this->mMaxMassFrc : nullptr);
    }

    const PhaseTransformationLaw::TKineticsType::Pointer PhaseTransformationLaw::GetKineticsModel()
    {
        return (this->mPhaseTranType == PhaseTransformationType::UNDEFINED_TRANSFORMATION) ? nullptr :this->mpKinetics;
    }

    const PhaseTransformationLaw::TPartitionType::Pointer PhaseTransformationLaw::GetPartitionModel()
    {
        return (this->mPhaseTranType == PhaseTransformationType::UNDEFINED_TRANSFORMATION || 
                this->mPhaseTranType == PhaseTransformationType::AUSTENIZATION || 
                this->mPhaseTranType == PhaseTransformationType::MARTENSITIC_TRANSFORMATION) ? nullptr :this->mpPartition;
    }

    const PhaseTransformationLaw::TTripType::Pointer PhaseTransformationLaw::GetTripModel()
    {
        return (this->mPhaseTranType == PhaseTransformationType::UNDEFINED_TRANSFORMATION) ? nullptr :this->mpTrip;
    }


    int PhaseTransformationLaw::Initialize(const Properties& rMaterialProperties, PhaseTransformationType rType)
    {
        int iret = 0x0;
        bool chk;

        this->mPhaseTranType = rType;

        KRATOS_TRY;

        if(this->mstat == Status::UnSet)
        {
            switch(mPhaseTranType)
            {
                case PhaseTransformationType::AUSTENIZATION:
                    // initialize the variable parameters
                    this->mLatent = &LATENT_AUSTENIZATION;
                    this->mTransStrain = &STRAIN_TRANSFORMATION_AUSTENIZATION;
                    this->mUpperT = &MAX_TEMPERATURE_AUSTENIZATION;
                    this->mLowerT = &MIN_TEMPERATURE_AUSTENIZATION;
                    this->mMaxMassFrc = &MAX_FRACTION_AUSTENIZATION;

                    // initialize the kinetics model
                    if(this->mpKinetics == nullptr)
                    {
                        this->mpKinetics = rMaterialProperties.GetValue(KINETICS_AUSTENIZATION);
                    }
                    this->mpKinetics->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // austenization has no partition model
                    if(this->mpPartition != nullptr)
                    {
                        this->mpPartition == nullptr;
                    }

                    // initialize the TRIP model
                    if(this->mpTrip == nullptr)
                    {
                        this->mpTrip = rMaterialProperties.GetValue(TRIP_AUSTENIZATION);
                    }
                    this->mpTrip->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // check all the pointers
                    chk = this->mLatent != nullptr && this->mTransStrain != nullptr && this->mUpperT != nullptr &&
                          this->mLowerT != nullptr && this->mMaxMassFrc != nullptr && this->mpKinetics != nullptr &&
                          this->mpTrip != nullptr;
                    break;
                case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                    this->mLatent = &LATENT_EUTECTOID_DECOMPOSITION;
                    this->mTransStrain = &STRAIN_TRANSFORMATION_EUTECTOID_DECOMPOSITION;
                    this->mUpperT = &MAX_TEMPERATURE_EUTECTOID_DECOMPOSITION;
                    this->mLowerT = &MIN_TEMPERATURE_EUTECTOID_DECOMPOSITION;
                    this->mMaxMassFrc = &MAX_FRACTION_EUTECTOID_DECOMPOSITION;

                    // initialize the kinetics model
                    if(this->mpKinetics == nullptr)
                    {
                        this->mpKinetics = rMaterialProperties.GetValue(KINETICS_EUTECTOID_DECOMPOSITION);
                    }
                    this->mpKinetics->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // initialize the partition model
                    if(this->mpPartition == nullptr)
                    {
                        this->mpPartition = rMaterialProperties.GetValue(PARTITION_EUTECTOID_DECOMPOSITION);
                    }
                    this->mpPartition->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // initialize the TRIP model
                    if(this->mpTrip == nullptr)
                    {
                        this->mpTrip = rMaterialProperties.GetValue(TRIP_EUTECTOID_DECOMPOSITION);
                    }
                    this->mpTrip->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // check all the pointers
                    chk = this->mLatent != nullptr && this->mTransStrain != nullptr && this->mUpperT != nullptr &&
                          this->mLowerT != nullptr && this->mMaxMassFrc != nullptr && this->mpKinetics != nullptr && 
                          this->mpPartition != nullptr && this->mpTrip != nullptr;
                    break;
                case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                    this->mLatent = &LATENT_EUTECTOID_DECOMPOSITION;
                    this->mTransStrain = &STRAIN_TRANSFORMATION_BAINITIC_TRANSFORMATION;
                    this->mUpperT = &MAX_TEMPERATURE_BAINITIC_TRANSFORMATION;
                    this->mLowerT = &MIN_TEMPERATURE_BAINITIC_TRANSFORMATION;
                    this->mMaxMassFrc = &MAX_FRACTION_BAINITIC_TRANSFORMATION;

                    // initialize the kinetics model
                    if(this->mpKinetics == nullptr)
                    {
                        this->mpKinetics = rMaterialProperties.GetValue(KINETICS_BAINITIC_TRANSFORMATION);
                    }
                    this->mpKinetics->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // initialize the partition model
                    if(this->mpPartition == nullptr)
                    {
                        this->mpPartition = rMaterialProperties.GetValue(PARTITION_BAINITIC_TRANSFORMATION);
                    }
                    this->mpPartition->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // initialize the TRIP model
                    if(this->mpTrip == nullptr)
                    {
                        this->mpTrip = rMaterialProperties.GetValue(TRIP_BAINITIC_TRANSFORMATION);
                    }
                    this->mpTrip->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // check all the pointers
                    chk = this->mLatent != nullptr && this->mTransStrain != nullptr && this->mUpperT != nullptr &&
                          this->mLowerT != nullptr && this->mMaxMassFrc != nullptr && this->mpKinetics != nullptr && 
                          this->mpPartition != nullptr && this->mpTrip != nullptr;
                    break;
                case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                    this->mLatent = &LATENT_EUTECTOID_DECOMPOSITION;
                    this->mTransStrain = &STRAIN_TRANSFORMATION_MARTENSITIC_TRANSFORMATION;
                    this->mUpperT = &MAX_TEMPERATURE_MARTENSITIC_TRANSFORMATION;
                    this->mLowerT = &MIN_TEMPERATURE_MARTENSITIC_TRANSFORMATION;
                    this->mMaxMassFrc = &MAX_FRACTION_MARTENSITIC_TRANSFORMATION;

                    // initialize the kinetics model
                    if(this->mpKinetics == nullptr)
                    {
                        this->mpKinetics = rMaterialProperties.GetValue(KINETICS_MARTENSITIC_TRANSFORMATION);
                    }
                    this->mpKinetics->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // martensitic transformation has no partition model
                    if(this->mpPartition != nullptr)
                    {
                        this->mpPartition = nullptr;
                    }

                    // initialize the TRIP model
                    if(this->mpTrip == nullptr)
                    {
                        this->mpTrip = rMaterialProperties.GetValue(TRIP_MARTENSITIC_TRANSFORMATION);
                    }
                    this->mpTrip->Initialize(rMaterialProperties, this->mPhaseTranType);

                    // check all the pointers
                    chk = this->mLatent != nullptr && this->mTransStrain != nullptr && this->mUpperT != nullptr &&
                          this->mLowerT != nullptr && this->mMaxMassFrc != nullptr && this->mpKinetics != nullptr && 
                          this->mpTrip != nullptr;
                    break;
                default:
                    KRATOS_ERROR<< "PhaseTransformationLaw:: The Phase Transformation Type is not defined." << std::endl;

            }
            if(chk)
            {
                this->mstat = Status::Set;
            }
            else
            {
                KRATOS_ERROR<< "Members for "<< this->GetPhaseTranTypeName()<<" are not assigned properly!" << std::endl;
            }

            iret = 0x1;

            KRATOS_INFO("PhaseTransformationLaw::")
                << this->GetPhaseTranTypeName() << " law is initialized. \n"
                << "  Kinetics model   : " << (this->GetKineticsModel() ? this->GetKineticsModel()->Info() : "None") << "\n"
                << "  Partition model  : " << (this->GetPartitionModel() ? this->GetPartitionModel()->Info() : "None") << "\n"
                << "  Trip model       : " << (this->GetTripModel() ? this->GetTripModel()->Info() : "None") << "\n"
                << std::endl;
        }

        KRATOS_CATCH("");
        return iret;
    }

    /**
     * @brief Returns the mean value integral value of unit latent of phase transformation
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rCurrentProcessInfo the current ProcessInfo instance
     * @return the mean value of the latent heat of phase transformation
     */
    double PhaseTransformationLaw::GetPhaseTransLatent(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId)
    {
        //返回相变潜热
        double latent_heat = 0.0;

        latent_heat = rMaterialProperties.GetValue(*mLatent, rElementGeometry[NodeId], rCurrentProcessInfo);
        return latent_heat;
    }


    double PhaseTransformationLaw::GetPhaseTransStrain(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      SizeType NodeId)
    {
        //返回相变应变

        double phase_trans_strain = 0.0;

        phase_trans_strain = rMaterialProperties.GetValue(*mTransStrain, rElementGeometry[NodeId], rCurrentProcessInfo);

        return phase_trans_strain;
    }

    double PhaseTransformationLaw::GetPhaseTransUpperT(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      SizeType NodeId)
    {
        //返回相变上限温度
        double upper_T = 0.0;

        upper_T = rMaterialProperties.GetValue(*mUpperT, rElementGeometry[NodeId], rCurrentProcessInfo);

        return upper_T;
    }

    double PhaseTransformationLaw::GetPhaseTransLowerT(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      SizeType NodeId)
    {
        //返回相变下限温度
        double lower_T = 0.0;

        lower_T = rMaterialProperties.GetValue(*mLowerT, rElementGeometry[NodeId], rCurrentProcessInfo);

        return lower_T;
    }

    double PhaseTransformationLaw::GetMaxTransMassFrc(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      SizeType NodeId)
    {
        //返回最大相变质量分数
        double max_frac = 0.0;
        max_frac = rMaterialProperties.GetValue(*mMaxMassFrc, rElementGeometry[NodeId], rCurrentProcessInfo);
        return max_frac;
    }

    int PhaseTransformationLaw::GetPartitionInfo(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double rMassFrac, double *rMassTransFracInc)
    {
        //返回相变配分信息，调用mPartition
        int iret = 0x0;
        if(this->mpPartition != nullptr)
        {
            iret = this->mpPartition->GetPartitionInfo(rMaterialProperties, rElementGeometry,
                                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId, rMassFrac, rMassTransFracInc);
        }
        return iret;
    }

    int PhaseTransformationLaw::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                        const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                        SizeType NodeId, double rIncuFrac, double rMassFrac, double rSatuFrac,
                                        double &rIncuFracInc, double &rMassFracInc)
    {
        //返回相变质量分数增量，调用mKinetics
        int iret = 0;

        if(this->mpKinetics != nullptr)
        {
             iret = this->mpKinetics->GetPhaseTransMassFrcInc(rMaterialProperties, rElementGeometry,
                                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId, rIncuFrac, rMassFrac, 
                                        rSatuFrac, rIncuFracInc, rMassFracInc);
        }
        return iret;
    }
    
    double PhaseTransformationLaw::GetEqTripInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                       const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                       SizeType NodeId, double rMassFrac, double rMassFracInc)
    {
        //返回等效Trip增量
        double dret = 0.0;

        if(this->mpTrip != nullptr)
        {
            dret = this->mpTrip->GetEqTripInc(rMaterialProperties, rElementGeometry,
                                        rShapeFunctionsValues, rCurrentProcessInfo, NodeId, rMassFrac, rMassFracInc);
        }
        else
        {
            KRATOS_ERROR << "mpTrip is not initialized in PhaseTransformationLaw::GetEqTripInc" << std::endl;
        }
        return dret;
    }

    int PhaseTransformationLaw::Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR <<  "Called the virtual function for Check from PhaseTransformationLaw"<< std::endl;
        return -1; // Return an error code
    }

    void PhaseTransformationLaw::GetModelFeatures(ModelFeatures& rModelFeatures)
    {
        KRATOS_ERROR <<  "Called the virtual function for GetModelFeatures from PhaseTransformationLaw"<< std::endl;
    }
}
