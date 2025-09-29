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

#include "kme.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Kme::Kme() 
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
    Kme::Kme(const Kme &rOther) 
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
    KineticsModel::Pointer Kme::Clone() const
    {
        return Kratos::make_shared<Kme>(Kme(*this));
    }

    /**
     * Destructor.
     */
    Kme::~Kme() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 2; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int Kme::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
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
                KRATOS_ERROR << "The Phase Transformation Type for Kme is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for Kme are not assigned" << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }

        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /** 
     * Kme 的参数1.critical_temperature、2.alpha
     */
    //Austenization kinetics
    int Kme::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
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
        double Tc, alpha, Tm, mass_frac_inc, time, dt;
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

                    Tc = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    alpha = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);

                    if (alpha > 0)
                    {
                        Tm = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                                    rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1));

                        if (Tm < Tc)
                        {
                            mass_frac_inc = 1 - exp(-alpha * (Tc - Tm)) - rMassFrac;

                            /* check the mass fraction increment */
                            if (mass_frac_inc > 0)
                            {
                                rMassFracInc = mass_frac_inc;
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
                                rIncuFrac = 1;
                            }
                        }
                        else
                        {
                            rMassFracInc = 0;
                            rIncuFracInc = 0;
                        }
                        iret = 0X1;
                    }
                }
            }
        }

        return iret; // Return 0 for successful calculation
    }
}
