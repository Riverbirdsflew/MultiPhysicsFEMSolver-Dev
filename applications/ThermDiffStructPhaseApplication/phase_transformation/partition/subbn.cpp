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

#include "subbn.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    SubBn::SubBn() 
        : PartitionModel()
    {
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;
        this->mPara = nullptr;
    }

    /**
     * Copy Constructor.
     */
    SubBn::SubBn(const SubBn &rOther) 
        : PartitionModel(rOther)
    {
        this->mType = rOther.mType;
        this->mstat = rOther.mstat;
        this->mPara = rOther.mPara;
    }

    /**
     * Clone Constructor.
     */
    PartitionModel::Pointer SubBn::Clone() const
    {
        return Kratos::make_shared<SubBn>(SubBn(*this));
    }

    /**
     * Destructor.
     */
    SubBn::~SubBn() 
    {
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;
        this->mPara = nullptr;
    }

    std::string SubBn::GetTransTypeName() const 
    {
        switch (mType) {
            case PhaseTransformationType::UNDEFINED_TRANSFORMATION:
                return "UNDEFINED_TRANSFORMATION";
            case PhaseTransformationType::AUSTENIZATION:
                return "AUSTENIZATION";
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                return "EUTECTOID_DECOMPOSITION";
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                return "BAINITIC_TRANSFORMATION";
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                return "MARTENSITIC_TRANSFORMATION";
            default:
                return "UNKNOWN_TRANSFORMATION";
        }
    }

    int SubBn::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
    {
        bool chk;
        int iret = 0x0;
        this->mType = rType;

        KRATOS_TRY;

        if (this->mstat == Status::UnSet)
        {
            switch (mType)
            {
            case PhaseTransformationType::AUSTENIZATION:
                KRATOS_ERROR << "The model " << this->Info() << " shoud not assigned to this type: " << this->GetTransTypeName() << std::endl;
                break;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                KRATOS_ERROR << "The model " << this->Info() << " shoud not assigned to this type: " << this->GetTransTypeName() << std::endl;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                this->mPara = &PARTITION_PARA_BAINITIC_TRANSFORMATION_1;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                KRATOS_ERROR << "The model " << this->Info() << " shoud not assigned to this type: " << this->GetTransTypeName() << std::endl;
                break;
            default:
                KRATOS_ERROR << "The Phase Transformation Type for  SubBn is not defined." << std::endl;
            }

            chk = this->mPara != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for SubBn are not assigned" << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }

        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /** 
     * 返回贝氏体分解配分结果：rMassTransFracInc = a/ferrite, 1 = cementite, 2 = pearlite, 3 = upper bainite, 4 = lower bainite
     * rMassFrac = partition mass fraction
     * Bainitic partition 的参数: 1.critical_temperature 
     */
    int SubBn::GetPartitionInfo(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                     const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                     SizeType NodeId, double rMassFrac, double *rMassTransFracInc)
    {
        
        int iret = 0X0;        
        double T_c, T_mean, dt, t;

        dt = rCurrentProcessInfo[DELTA_TIME];
        t = rCurrentProcessInfo[TIME];

        if (t >= 0 && dt >= 0 && rMassFrac > 0 && rMassTransFracInc != nullptr)
        {
            for (SizeType i = 0; i < 2; i++)
            {
                rMassTransFracInc[i] = 0;
            }

            T_c = rMaterialProperties.GetValue(*mPara, rElementGeometry[NodeId], rCurrentProcessInfo);
            T_mean = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                            rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1));

            if (T_mean > T_c)
            {
                rMassTransFracInc[0] = rMassFrac;
            }
            else
            {
                rMassTransFracInc[1] = rMassFrac;
            }
            iret = 0x1;
        }

        return iret; // Return 0 for successful calculation
    }
}
