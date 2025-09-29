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

#include "maxwell_power.h"
#include "therm_diff_struct_phase_application_variables.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    MWPower::MWPower() 
    {
        unsigned int i;

        this->mType = PhaseType::UNDEFINED;
        this->mstat = Status::UnSet;

        for (i = 0; i < 5; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    /**
     * Copy Constructor.
     */
    MWPower::MWPower(const MWPower &rOther) 
    {
        unsigned int i;

        this->mstat = rOther.mstat;
        this->mType = rOther.mType;
        
        for (i = 0; i < 5; i++)
        {
            this->mPara[i] = rOther.mPara[i];
        }
    }

    /**
     * Clone Constructor.
     */
    CreepModel::Pointer MWPower::Clone() const
    {
        return Kratos::make_shared<MWPower>(MWPower(*this));
    }

    /**
     * Destructor.
     */
    MWPower::~MWPower() 
    {
        unsigned int i;

        this->mstat = Status::UnSet;
        this->mType = PhaseType::UNDEFINED;
        
        for (i = 0; i < 5; i++)
        {
            this->mPara[i] = nullptr;
        }
    }

    int MWPower::Initialize(const Properties &rMaterialProperties, PhaseType rType)
    {
        bool chk;
        int iret = 0x0;

        this->mType = rType;

        KRATOS_TRY;

        if (this->mstat == Status::UnSet)
        {
            switch (mType)
            {
            case PhaseType::AUSTENITE:
                this->mPara[0] = &CREEP_PARA_AUSTENITE_1;
                this->mPara[1] = &CREEP_PARA_AUSTENITE_2;
                this->mPara[2] = &CREEP_PARA_AUSTENITE_3;
                this->mPara[1] = &CREEP_PARA_AUSTENITE_4;
                this->mPara[2] = &CREEP_PARA_AUSTENITE_5;
                break;
            case PhaseType::FERRITE:
                this->mPara[0] = &CREEP_PARA_FERRITE_1;
                this->mPara[1] = &CREEP_PARA_FERRITE_2;
                this->mPara[2] = &CREEP_PARA_FERRITE_3;
                this->mPara[3] = &CREEP_PARA_FERRITE_4;
                this->mPara[4] = &CREEP_PARA_FERRITE_5;
                break;
            case PhaseType::CEMENTITE:
                this->mPara[0] = &CREEP_PARA_CEMENTITE_1;
                this->mPara[1] = &CREEP_PARA_CEMENTITE_2;
                this->mPara[2] = &CREEP_PARA_CEMENTITE_3;
                this->mPara[3] = &CREEP_PARA_CEMENTITE_4;
                this->mPara[4] = &CREEP_PARA_CEMENTITE_5;
                break;
            case PhaseType::PEARLITE:
                this->mPara[0] = &CREEP_PARA_PEARLITE_1;
                this->mPara[1] = &CREEP_PARA_PEARLITE_2;
                this->mPara[2] = &CREEP_PARA_PEARLITE_3;
                this->mPara[3] = &CREEP_PARA_PEARLITE_4;
                this->mPara[4] = &CREEP_PARA_PEARLITE_5;
                break;
            case PhaseType::UPPER_BAINITE:
                this->mPara[0] = &CREEP_PARA_UPPER_BAINITE_1;
                this->mPara[1] = &CREEP_PARA_UPPER_BAINITE_2;
                this->mPara[2] = &CREEP_PARA_UPPER_BAINITE_3;
                this->mPara[3] = &CREEP_PARA_UPPER_BAINITE_4;
                this->mPara[4] = &CREEP_PARA_UPPER_BAINITE_5;
                break;
            case PhaseType::LOWER_BAINITE:
                this->mPara[0] = &CREEP_PARA_LOWER_BAINITE_1;
                this->mPara[1] = &CREEP_PARA_LOWER_BAINITE_2;
                this->mPara[2] = &CREEP_PARA_LOWER_BAINITE_3;
                this->mPara[3] = &CREEP_PARA_LOWER_BAINITE_4;
                this->mPara[4] = &CREEP_PARA_LOWER_BAINITE_5;
                break;
            case PhaseType::MARTENSITE:
                this->mPara[0] = &CREEP_PARA_MARTENSITE_1;
                this->mPara[1] = &CREEP_PARA_MARTENSITE_2;
                this->mPara[2] = &CREEP_PARA_MARTENSITE_3;
                this->mPara[3] = &CREEP_PARA_MARTENSITE_4;
                this->mPara[4] = &CREEP_PARA_MARTENSITE_5;
            default:
                KRATOS_ERROR << "The Phase Type for Maxwell Linear is not defined." << std::endl;
            }

            chk = this->mPara[0] != nullptr && this->mPara[1] != nullptr && this->mPara[2] != nullptr && this->mPara[3] != nullptr && this->mPara[4] != nullptr;
            if (!chk)
                KRATOS_ERROR << "The parameters for Maxwell Linear are not assigned." << std::endl;
            else
                this->mstat = Status::Set;

            iret = 0x1;
        }
        KRATOS_CATCH("");
        return iret; // Return 0 for successful initialization
    }

    /**
     * 返回指数蠕变模型的参数个数和数组
     * creep的参数1.factor_a 2.factor_stress、3.factor_creep 4.factor_temperature  5.factor_time
     */
    // creep parameters
    int MWPower::GetCreepParameters(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                   const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo,
                                   SizeType NodeId, SizeType *pnm, double *pdata)
    {
        int iret = 0x0;
        unsigned int i;
        double cp[5], dt, t;

        dt = rCurrentProcessInfo[DELTA_TIME];
        t = rCurrentProcessInfo[TIME];

        if (NodeId > 0 && dt >= 0 && t >= 0 && pnm != nullptr && pdata != nullptr)
        {
            if (this->mstat == Status::Set)
            {
                cp[0] = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[1] = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[2] = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[3] = rMaterialProperties.GetValue(*(mPara[3]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[4] = rMaterialProperties.GetValue(*(mPara[4]), rElementGeometry[NodeId], rCurrentProcessInfo);

                for (i = 0; i < 5; i++)
                {
                    pdata[i] = cp[i];
                }

                *pnm = 5; // Number of parameters
                iret = 0x1;
            }
        }
        return iret;
    }

    /**
     * 返回蠕变模型的等效参数
     * creep的参数1.factor_a 2.factor_stress、3.factor_creep 4.factor_temperature  5.factor_time
     */
    double MWPower::GetEqCreepRate(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                   const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                   SizeType NodeId)
    {
        unsigned int i;
        double dret, t;
        double cp[5];

        t = rCurrentProcessInfo[TIME];

        if (NodeId > 0 && t >= 0)
        {
            if (this->mstat == Status::Set)
            {
                cp[0] = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[1] = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[2] = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[3] = rMaterialProperties.GetValue(*(mPara[3]), rElementGeometry[NodeId], rCurrentProcessInfo);
                cp[4] = rMaterialProperties.GetValue(*(mPara[4]), rElementGeometry[NodeId], rCurrentProcessInfo);

                double eq_s = rElementGeometry[NodeId].FastGetSolutionStepValue(EQUIVALENT_STRESS);
                double eq_creep_e = rElementGeometry[NodeId].FastGetSolutionStepValue(EQUIVALENT_STRAIN_CREEP);
                double Tem = rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE);
                cp[1] = pow(eq_s, cp[1]);
                cp[2] = pow(eq_creep_e, cp[2]);
                cp[3] = pow(Tem, cp[3]);
                cp[4] = pow(t, cp[4] - 1);

                dret = cp[0] * cp[1] * cp[2] * cp[3] * cp[4];
            }
        }

        return dret;
    }

    double MWPower::GetEqCreepRateMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                       const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                       SizeType NodeId)
    {
        unsigned int i;
        double dret, t, dt;
        double cp1[5], cp2[5];
        t = rCurrentProcessInfo[TIME];
        dt = rCurrentProcessInfo[DELTA_TIME];

        if (NodeId <= 0 && t >= 0 && dt >= 0)
        {
            if (this->mstat == Status::Set)
            {
                cp1[0] = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], 0, rCurrentProcessInfo);
                cp1[1] = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], 0, rCurrentProcessInfo);
                cp1[2] = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], 0, rCurrentProcessInfo);
                cp1[3] = rMaterialProperties.GetValue(*(mPara[3]), rElementGeometry[NodeId], 0, rCurrentProcessInfo);
                cp1[4] = rMaterialProperties.GetValue(*(mPara[4]), rElementGeometry[NodeId], 0, rCurrentProcessInfo);
                cp2[0] = rMaterialProperties.GetValue(*(mPara[0]), rElementGeometry[NodeId], 1, rCurrentProcessInfo);
                cp2[1] = rMaterialProperties.GetValue(*(mPara[1]), rElementGeometry[NodeId], 1, rCurrentProcessInfo);
                cp2[2] = rMaterialProperties.GetValue(*(mPara[2]), rElementGeometry[NodeId], 1, rCurrentProcessInfo);
                cp2[3] = rMaterialProperties.GetValue(*(mPara[3]), rElementGeometry[NodeId], 1, rCurrentProcessInfo);
                cp2[4] = rMaterialProperties.GetValue(*(mPara[4]), rElementGeometry[NodeId], 1, rCurrentProcessInfo);

                for (i = 0; i < 5; i++)
                {
                    cp2[i] = (cp1[i] + cp2[i]) / 2.0; // Mean value
                }
                double eq_s = rElementGeometry[NodeId].FastGetSolutionStepValue(EQUIVALENT_STRESS);
                double eq_creep_e = rElementGeometry[NodeId].FastGetSolutionStepValue(EQUIVALENT_STRAIN_CREEP);
                double Tem = rElementGeometry[NodeId].FastGetSolutionStepValue(TEMPERATURE);
                cp2[1] = pow(eq_s, cp2[1]);
                cp2[2] = pow(eq_creep_e, cp2[2]);
                cp2[3] = pow(Tem, cp2[3]);
                cp2[4] = pow(t, cp2[4] - 1);

                dret = cp2[0] * cp2[1] * cp2[2] * cp2[3] * cp2[4];
            }

            return dret;
        }
    }
}
