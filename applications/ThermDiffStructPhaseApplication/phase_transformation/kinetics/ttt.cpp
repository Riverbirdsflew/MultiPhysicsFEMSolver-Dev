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

#include "ttt.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Ttt::Ttt() 
    {
        unsigned int i;

        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;

        for (i = 0; i < 2; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    /**
     * Copy Constructor.
     */
    Ttt::Ttt(const Ttt &rOther) 
         : mType(rOther.mType)
    {
        unsigned int i;

        this->mstat = rOther.mstat;
        
        for (i = 0; i < 2; i++)
        {
            this->mPara[i] = rOther.mPara[i];
        }
    }

    /**
     * Clone Constructor.
     */
    KineticsModel::Pointer Ttt::Clone() const
    {
        return Kratos::make_shared<Ttt>(Ttt(*this));
    }

    /**
     * Destructor.
     */
    Ttt::~Ttt() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 2; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int Ttt::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
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
                break;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                this->mPara[0] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_1;
                this->mPara[1] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_2;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_2;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_2;
            default:
                KRATOS_ERROR << "The Phase Transformation Type for TTT is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr && this->mPara[2] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for TTT are not assigned" << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }
        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /** 
     * TTT 的参数1.incubation、2.completed_time
     */
    int Ttt::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
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
        double T1, T99, Tm, n, b, t_tf, dt, time;
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

                    T1 = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    T99 = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);

                    if (T1 > 0 && T99 > T1)
                    {
                        n = PAR_CONSTANT_CN / log10(T99 / T1);
                        b = PAR_CONSTANT_CB / pow(T99, n);

                        t_tf = pow(-log(1 - rMassFrac) / b, 1 / n);
                        rMassFracInc = 1 - exp(-b * pow(t_tf + dt, n)) - rMassFrac;

                        /* get the cfinc */
                        if (rIncuFrac >= 1)
                        {
                            rIncuFracInc = 0;
                        }
                        else
                        {
                            rIncuFracInc = 1;
                        }
                        iret = 0X1;
                    }
                }
            }
        }

        return iret; // Return 0 for successful calculation
    }
}
