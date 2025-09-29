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

#include "iso.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Iso::Iso()
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
    Iso::Iso(const Iso &rOther) 
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
    KineticsModel::Pointer Iso::Clone() const
    {
        return Kratos::make_shared<Iso>(Iso(*this));
    }

    /**
     * Destructor.
     */
    Iso::~Iso() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 2; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int Iso::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
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
                KRATOS_ERROR << "The Phase Transformation Type for Iso is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for Iso are not assigned" << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }

        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /**
     * Iso 的参数1.factor_a、2.factor_e
     */
    //Austenization kinetics
    int Iso::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
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
        double A, E, T_mean, rate, remain, t, dt;
        dt = rCurrentProcessInfo[DELTA_TIME];
        t = rCurrentProcessInfo[TIME];

        A = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
        E = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
        T_mean = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                        rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1)) +
                 273.15;
        if (NodeId > 0)
        {
            if (dt >= 0 && time >= 0 && rIncuFrac >= 0 && rMassFrac >= 0)
            {
                if (this->mstat == Status::Set)
                {
                    if (A > 0 && E > 0)
                    {
                        rate = A * exp(-E / (PAR_CONSTANT_R * T_mean)) * dt;
                        remain = 1.0 - rMassFrac;

                        if (rate > remain)
                        {
                            rMassFracInc = remain;
                        }
                        else
                        {
                            rMassFracInc = rate;
                        }

                        // incubation fraction
                        if (rIncuFrac >= 1)
                        {
                            rIncuFracInc = 0.0;
                        }
                        else
                        {
                            rIncuFracInc = 1.0;
                        }
                        iret = 0X1;
                    }
                }
            }
        }

        return iret; // Return 0 for successful calculation
    }
}
