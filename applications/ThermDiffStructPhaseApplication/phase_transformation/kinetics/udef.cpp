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

#include "udef.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Udef::Udef()
    {
        unsigned int i;

        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;

        for (i = 0; i < 11; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    /**
     * Copy Constructor.
     */
    Udef::Udef(const Udef &rOther)
        : mType(rOther.mType)
    {
        unsigned int i;

        this->mstat = rOther.mstat;
        
        for (i = 0; i < 11; i++)
        {
            this->mPara[i] = rOther.mPara[i];
        }
    }

    /**
     * Clone Constructor.
     */
    KineticsModel::Pointer Udef::Clone() const
    {
        return Kratos::make_shared<Udef>(Udef(*this));
    }

    /**
     * Destructor.
     */
    Udef::~Udef() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 11; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int Udef::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
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
                this->mPara[0] = &KINETICS_PARA_AUSTENIZATION_1;
                this->mPara[1] = &KINETICS_PARA_AUSTENIZATION_2;
                this->mPara[2] = &KINETICS_PARA_AUSTENIZATION_3;
                this->mPara[3] = &KINETICS_PARA_AUSTENIZATION_4;
                this->mPara[4] = &KINETICS_PARA_AUSTENIZATION_5;
                this->mPara[5] = &KINETICS_PARA_AUSTENIZATION_6;
                this->mPara[6] = &KINETICS_PARA_AUSTENIZATION_7;
                this->mPara[7] = &KINETICS_PARA_AUSTENIZATION_8;
                this->mPara[8] = &KINETICS_PARA_AUSTENIZATION_9;
                this->mPara[9] = &KINETICS_PARA_AUSTENIZATION_10;
                this->mPara[10] = &KINETICS_PARA_AUSTENIZATION_11;
                break;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                this->mPara[0] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_1;
                this->mPara[1] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_2;
                this->mPara[2] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_3;
                this->mPara[3] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_4;
                this->mPara[4] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_5;
                this->mPara[5] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_6;
                this->mPara[6] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_7;
                this->mPara[7] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_8;
                this->mPara[8] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_9;
                this->mPara[9] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_10;
                this->mPara[10] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_11;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_2;
                this->mPara[2] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_3;
                this->mPara[3] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_4;
                this->mPara[4] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_5;
                this->mPara[5] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_6;
                this->mPara[6] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_7;
                this->mPara[7] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_8;
                this->mPara[8] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_9;
                this->mPara[9] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_10;
                this->mPara[10] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_11;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_2;
                this->mPara[2] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_3;
                this->mPara[3] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_4;
                this->mPara[4] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_5;
                this->mPara[5] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_6;
                this->mPara[6] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_7;
                this->mPara[7] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_8;
                this->mPara[8] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_9;
                this->mPara[9] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_10;
                this->mPara[10] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_11;
            default:
                KRATOS_ERROR << "The Phase Transformation Type for Udef is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr && this->mPara[2] != nullptr && this->mPara[3] != nullptr && this->mPara[4] != nullptr && this->mPara[5] != nullptr && this->mPara[6] != nullptr && this->mPara[7] != nullptr && this->mPara[8] != nullptr && this->mPara[9] != nullptr && this->mPara[10] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for Udef are not assigned." << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }

        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /** 
     * Udef 的参数1.incubation、2.10%_time 3.20%_time 4.30%_time 5.40%_time 6.50%_time 7.60%_time 8.70%_time 9.80%_time 10.90%_time 11.complete_time
     */
    //Austenization kinetics
    int Udef::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
                                      const GeometryType &rElementGeometry,
                                      const Vector &rShapeFunctionsValues,
                                      const ProcessInfo &rCurrentProcessInfo,
                                      SizeType NodeId,
                                      double rIncuFrac,
                                      double rMassFrac,
                                      double rSatuFrac,
                                      double &rIncuFracInc,
                                      double &rMassFracInc)
    {
        int iret = 0X0;
        double left, right, inter_t, inter_rMassFrac, dt, time;
        double incu, para[11];
        int i, j;
        dt = rCurrentProcessInfo[DELTA_TIME];
        time = rCurrentProcessInfo[TIME];

        if (NodeId > 0)
        {
            if (dt >= 0 && time >= 0 && rIncuFrac >= 0 && rMassFrac >= 0)
            {
                if (this->mstat == Status::Set)
                {
                    rIncuFracInc = 0.0;
                    rMassFracInc = 0.0;

                    para[0] = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[1] = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[2] = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[3] = rMaterialProperties.GetValue(*(mPara[3]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[4] = rMaterialProperties.GetValue(*(mPara[4]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[5] = rMaterialProperties.GetValue(*(mPara[5]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[6] = rMaterialProperties.GetValue(*(mPara[6]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[7] = rMaterialProperties.GetValue(*(mPara[7]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[8] = rMaterialProperties.GetValue(*(mPara[8]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[9] = rMaterialProperties.GetValue(*(mPara[9]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    para[10] = rMaterialProperties.GetValue(*(mPara[10]), rElementGeometry[NodeId], rCurrentProcessInfo);

                    /* check the ts and cf */
                    if (para[0] <= 0 || rIncuFrac >= 1)
                    {
                        /* search the rMassFrac */
                        for (i = 0, left = 0; i < 10; i++)
                        {
                            /* get tne right limit */
                            right = left + 0.1;

                            /* check the rMassFrac */
                            if (rMassFrac >= left && rMassFrac <= right)
                            {
                                /* get the virtual time */
                                inter_t = 10 * ((rMassFrac - left) * para[i + 1] + (right - rMassFrac) * para[i]) + dt;

                                /* exit cycle */
                                break;
                            }

                            left = right;
                        }

                        /* calculate the mass fraction increment directly */
                        for (; i < 10; i++)
                        {
                            j = i + 1;

                            /* check the vritual time */
                            if (inter_t >= para[i] && inter_t <= para[j])
                            {
                                /* get the mass fraction */
                                inter_rMassFrac = 0.1 * (i + (inter_t - para[i]) / (para[j] - para[i]));

                                /* exit cycle */
                                break;
                            }
                        }

                        /* check the i */
                        if (i == 10)
                        {
                            inter_rMassFrac = 1;
                        }

                        /* check the increment */
                        if (inter_rMassFrac > rMassFrac)
                        {
                            /* get the mass fraction increment */
                            rMassFracInc = inter_rMassFrac - rMassFrac;
                        }
                        else
                        {
                            rMassFracInc = 0;
                        }

                        /* get the cfinc */
                        if (rIncuFrac >= 1)
                        {
                            rIncuFracInc = 0;
                        }
                        else
                        {
                            rIncuFracInc = 1;
                        }
                    }
                    else
                    {
                        /* calculate the incubation fraction increment */
                        inter_rMassFrac = dt / para[0];

                        /* check the incubation fraction */
                        if (inter_rMassFrac + rIncuFrac > 1)
                        {
                            /* calculate the effective time increment */
                            inter_rMassFrac = para[0] * (inter_rMassFrac + rIncuFrac - 1);

                            /* search the rMassFrac */
                            for (i = 0, left = 0; i < 10; i++)
                            {
                                /* get tne right limit */
                                right = left + 0.1;

                                /* check the rMassFrac */
                                if (rMassFrac >= left && rMassFrac <= right)
                                {
                                    /* get the virtual time */
                                    inter_t = 10 * ((rMassFrac - left) * para[i + 1] + (right - rMassFrac) * para[i]) + inter_rMassFrac;

                                    /* exit cycle */
                                    break;
                                }
                            }

                            /* calculate the mass fraction increment directly */
                            for (; i < 10; i++)
                            {
                                j = i + 1;

                                /* check the vritual time */
                                if (inter_t >= para[i] && inter_t <= para[j])
                                {
                                    /* get the mass fraction */
                                    inter_rMassFrac = 0.1 * (i + (inter_t - para[i]) / (para[j] - para[i]));

                                    /* exit cycle */
                                    break;
                                }
                            }

                            /* check the i */
                            if (i == 10)
                            {
                                inter_rMassFrac = 1;
                            }

                            /* check the increment */
                            if (inter_rMassFrac > rMassFrac)
                            {
                                /* get the mass fraction increment */
                                rMassFracInc = inter_rMassFrac - rMassFrac;
                            }
                            else
                            {
                                rMassFracInc = 0;
                            }

                            /* cfinc = 1 - cf */
                            rIncuFracInc = 1 - rIncuFrac;
                        }
                        else
                        {
                            /* rMassFracinc = 0, cfinc = timinc / ts */
                            rMassFracInc = 0;
                            rIncuFracInc = inter_rMassFrac;
                        }
                    }
                    iret = 0X1;
                }
            }
        }

        return iret; // Return 0 for successful calculation
    }
}
