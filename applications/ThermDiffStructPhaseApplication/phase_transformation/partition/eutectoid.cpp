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

#include "eutectoid.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Eutectoid::Eutectoid() 
        : PartitionModel()
    {
        unsigned int i;

        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;

        for (i = 0; i < 3; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    /**
     * Copy Constructor.
     */
    Eutectoid::Eutectoid(const Eutectoid &rOther) 
        : mType(rOther.mType)
    {
        unsigned int i;

        this->mstat = rOther.mstat;
        
        for (i = 0; i < 3; i++)
        {
            this->mPara[i] = rOther.mPara[i];
        }
    }

    /**
     * Clone Constructor.
     */
    PartitionModel::Pointer Eutectoid::Clone() const
    {
        return Kratos::make_shared<Eutectoid>(Eutectoid(*this));
    }

    /**
     * Destructor.
     */
    Eutectoid::~Eutectoid() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 3; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    std::string Eutectoid::GetTransTypeName() const 
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


    int Eutectoid::Initialize(const Properties& rMaterialProperties, PhaseTransformationType rType)
    {
        bool chk;
        int iret = 0x0;
        
        this->mType = rType;

        KRATOS_TRY;

        if (this->mstat == Status::UnSet)
        {
            switch(mType)
        {
            case PhaseTransformationType::AUSTENIZATION:
                KRATOS_ERROR <<  "The model " << this->Info()<< " shoud not assigned to this type: " << this->GetTransTypeName() << std::endl;
                break;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                this->mPara[0] = &PARTITION_PARA_EUTECTOID_DECOMPOSITION_1;
                this->mPara[1] = &PARTITION_PARA_EUTECTOID_DECOMPOSITION_2;
                this->mPara[2] = &PARTITION_PARA_EUTECTOID_DECOMPOSITION_3;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                KRATOS_ERROR <<  "The model " << this->Info()<< " shoud not assigned to this type: " << this->GetTransTypeName() << std::endl;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                KRATOS_ERROR <<  "The model " << this->Info()<< " shoud not assigned to this type: " << this->GetTransTypeName() << std::endl;
                break;
            default:
                KRATOS_ERROR <<  "The Phase Transformation Type for Eutectoid is not defined."<< std::endl;
        }

        chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr && this->mPara[2] != nullptr;
        if (!chk)
            KRATOS_ERROR <<  "The parameters for Eutectoid are not assigned"<< std::endl;
        else
            this->mstat = Status::Set;

        iret = 0x1;
        }
        KRATOS_CATCH("");
        return iret; // Return 0 for successful initialization
    }

    /** 
     * 返回共析分解配分结果：rMassTransFracInc = a/ferrite, 1 = cementite, 2 = pearlite, 3 = upper bainite, 4 = lower bainite
     * rMassFrac = austenite mass fraction
     * Eutetoid 的参数: 1.eutectoid_temperature 2. min_saturated_carboncontent 3. max_saturated_carboncontant
     */
    int Eutectoid::GetPartitionInfo(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                    const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                    SizeType NodeId, double rMassFrac, double *rMassTransFracInc)
    {

        int iret = 0X0;
        double T_eute, min_c, max_c, c_alpha, T_mean, t, dt;

        dt = rCurrentProcessInfo[DELTA_TIME];
        t = rCurrentProcessInfo[TIME];

        if (NodeId > 0)
        {
            if (t >= 0 && dt >= 0 && rMassFrac > 0 && rMassTransFracInc != nullptr)
            {
                if (this->mstat != Status::Set)
                {
                    /* init the mtfinc */
                    for (SizeType i = 0; i < 3; i++)
                    {
                        rMassTransFracInc[i] = 0;
                    }

                    T_eute = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    min_c = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    max_c = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], rCurrentProcessInfo);

                    if (min_c > 0 && min_c <= max_c && max_c < 6.69)
                    {
                        c_alpha = rElementGeometry[NodeId].FastGetSolutionStepValue(CARBON_CONTENT_AUSTENITE);
                        T_mean = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                                        rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1));

                        if (T_mean > T_eute)
                        {
                            if (c_alpha >= 0 && c_alpha < min_c)
                            {
                                rMassTransFracInc[0] = (1 - c_alpha / min_c) * rMassFrac;
                            }
                            if (c_alpha >= min_c && c_alpha <= max_c)
                            {
                                rMassTransFracInc[1] = (c_alpha - max_c) / (6.69 - max_c) * rMassFrac;
                            }
                        }
                        else
                        {
                            if (c_alpha >= 0 && c_alpha < min_c)
                            {
                                rMassTransFracInc[0] = (1 - c_alpha / min_c) * rMassFrac;
                                rMassTransFracInc[2] = rMassFrac - rMassTransFracInc[0];
                            }
                            if (c_alpha >= min_c && c_alpha <= max_c)
                            {
                                rMassTransFracInc[1] = (c_alpha - max_c) / (6.69 - max_c) * rMassFrac;
                                rMassTransFracInc[2] = rMassFrac - rMassTransFracInc[1];
                            }
                        }
                        iret = 0X1;
                    }
                }
            }
        }

        return iret; // Return 0 for successful calculation
    }
}
