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

#include "avrami.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Avrami::Avrami()
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
    Avrami::Avrami(const Avrami &rOther)
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
    KineticsModel::Pointer Avrami::Clone() const
    {
        return Kratos::make_shared<Avrami>(Avrami(*this));
    }

    /**
     * Destructor.
     */
    Avrami::~Avrami() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 3; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int Avrami::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
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
                break;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                this->mPara[0] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_1;
                this->mPara[1] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_2;
                this->mPara[2] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_3;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_2;
                this->mPara[2] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_3;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_2;
                this->mPara[2] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_3;
            default:
                KRATOS_ERROR << "The Phase Transformation Type for Avrami is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr && this->mPara[2] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for Avrami are not assigned" << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }

        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /**
     * Avrami的参数1.incubation、2.factor_b、3.factor_n
     */
    //Austenization kinetics
    int Avrami::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
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
        double ts, b, n, t, t_rIncuFrac, dt, time;
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
                    
                    ts = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    b = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    n = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], rCurrentProcessInfo);

                    if (b > 0 && n > 0)
                    {
                        if (ts <= 0.0 || rIncuFrac >= 1.0)
                        {
                            // No incubation or already completed
                            t = pow(-log(1.0 - rMassFrac) / b, 1.0 / n);
                            rMassFracInc = 1 - exp(-b * pow(t + dt, n)) - rMassFrac;

                            if (rIncuFrac >= 1)
                            {
                                // Incubation complete
                                rIncuFracInc = 0.0;
                            }
                            else
                            {
                                rIncuFracInc = 1.0;
                            }
                        }
                        else
                        {
                            // Still in incubation
                            t_rIncuFrac = ts * (1.0 - rIncuFrac);

                            if (dt > t_rIncuFrac)
                            {
                                // If the time step exceeds the incubation time
                                t = pow(-log(1.0 - rMassFrac) / b, 1.0 / n);
                                rMassFracInc = 1 - exp(-b * pow(t + t_rIncuFrac, n)) - rMassFrac;
                                rIncuFracInc = 1.0 - rIncuFrac; // Incubation complete
                            }
                            else
                            {
                                // If the time step is within the incubation time
                                rMassFracInc = 0.0;
                                rIncuFracInc = dt / ts; // Still in incubation
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
