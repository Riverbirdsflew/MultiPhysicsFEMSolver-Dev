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

#include "jmak.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Jmak::Jmak() 
    {
        unsigned int i;

        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;

        for (i = 0; i < 5; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    /**
     * Copy Constructor.
     */
    Jmak::Jmak(const Jmak &rOther) 
        : mType(rOther.mType)
    {
        unsigned int i;

        this->mstat = rOther.mstat;
        
        for (i = 0; i < 5; i++)
        {
            this->mPara[i] = rOther.mPara[i];
        }
    }

    /**
     * Clone Constructor.
     */
    KineticsModel::Pointer Jmak::Clone() const
    {
        return Kratos::make_shared<Jmak>(Jmak(*this));
    }

    /**
     * Destructor.
     */
    Jmak::~Jmak() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        
        for (i = 0; i < 5; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int Jmak::Initialize(const Properties &rMaterialProperties, PhaseTransformationType rType)
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
                break;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                this->mPara[0] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_1;
                this->mPara[1] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_2;
                this->mPara[2] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_3;
                this->mPara[3] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_4;
                this->mPara[4] = &KINETICS_PARA_EUTECTOID_DECOMPOSITION_5;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_2;
                this->mPara[2] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_3;
                this->mPara[3] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_4;
                this->mPara[4] = &KINETICS_PARA_BAINITIC_TRANSFORMATION_5;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                this->mPara[0] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_1;
                this->mPara[1] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_2;
                this->mPara[2] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_3;
                this->mPara[3] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_4;
                this->mPara[4] = &KINETICS_PARA_MARTENSITIC_TRANSFORMATION_5;
                break;
            default:
                KRATOS_ERROR << "The Phase Transformation Type for Jmak is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr && this->mPara[2] != nullptr && this->mPara[3] != nullptr && this->mPara[4] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for Jmak are not assigned" << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }
        KRATOS_CATCH("");

        return iret; // Return 0 for successful initialization
    }

    /**
     * Jmak 的参数1.factor_kd、2.factor_q 、3.critical_temperature、4.factor_m、5.factor_n
     */
    //Austenization kinetics
    int Jmak::GetPhaseTransMassFrcInc(const Properties &rMaterialProperties,
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
        double K, Q, Tc, m, n, Tm, A, t, dt;

        if (NodeId > 0)
        {
            if (dt >= 0 && time >= 0 && rIncuFrac >= 0 && rMassFrac >= 0)
            {
                if (this->mstat == Status::Set)
                {
                    rIncuFracInc = 0.0;
                    rMassFracInc = 0.0;

                    K = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    Q = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    Tc = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    m = rMaterialProperties.GetValue(*(mPara[3]), rElementGeometry[NodeId], rCurrentProcessInfo);
                    n = rMaterialProperties.GetValue(*(mPara[4]), rElementGeometry[NodeId], rCurrentProcessInfo);

                    if (K > 0 && Q >= 0 && m >= 0 && n > 0)
                    {
                        Tm = 0.5 * (rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE) +
                                    rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE, 1)) +273.15;

                        if (Tm < Tc)
                        {
                            A = K * (Tc - Tm) * exp(-Q / (PAR_CONSTANT_R * Tm));
                            t = pow(-log(1 - rMassFrac) / A, 1 / n);
                            rMassFracInc = 1 - exp(-A * pow(t + dt, n)) - rMassFrac;
                        }
                        else
                        {
                            rMassFracInc = 0.0;
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
