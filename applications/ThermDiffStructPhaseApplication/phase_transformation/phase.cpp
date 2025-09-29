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

#include "phase.h"
#include "therm_diff_struct_phase_application_variables.h"
#include "includes/properties.h"

namespace Kratos
{
    /**
     * Constructor.
     */
    Phase::Phase() 
    {
        this->mstat = Status::UnSet;
        this->mPhaseType = PhaseType::UNDEFINED;
        this->mpDensity = nullptr;
        this->mpConductivity = nullptr;
        this->mpSpecificHeat = nullptr;
        this->mpDiffCond = nullptr;
        this->mpEModulus = nullptr;
        this->mpPoisson = nullptr;
        this->mpThermoExpansion = nullptr;
        this->mpYieldStrength = nullptr;
        this->mpCreep = nullptr;

    }

    Phase::Phase(PhaseType type, TCreepType::Pointer pCreepModel)
    {
        this->mstat = Status::UnSet;
        this->mPhaseType = type;
        this->mpCreep = pCreepModel;
        this->mpDensity = nullptr;
        this->mpConductivity = nullptr;
        this->mpSpecificHeat = nullptr;
        this->mpDiffCond = nullptr;
        this->mpEModulus = nullptr;
        this->mpPoisson = nullptr;
        this->mpThermoExpansion = nullptr;
        this->mpYieldStrength = nullptr;
    }

    Phase::Phase(const Phase &other)
    {
        this->mstat = other.mstat;
        this->mPhaseType = other.mPhaseType;
        this->mpCreep = other.mpCreep;
        this->mpDensity = other.mpDensity;
        this->mpConductivity = other.mpConductivity;
        this->mpSpecificHeat = other.mpSpecificHeat;
        this->mpDiffCond = other.mpDiffCond;
        this->mpEModulus = other.mpEModulus;
        this->mpPoisson = other.mpPoisson;
        this->mpThermoExpansion = other.mpThermoExpansion;
        this->mpYieldStrength = other.mpYieldStrength;
    }

    Phase::~Phase()
    {
        this->mstat = Status::UnSet;
        this->mpCreep = nullptr;
        this->mPhaseType = PhaseType::UNDEFINED;
        this->mpDensity = nullptr;
        this->mpConductivity = nullptr;
        this->mpSpecificHeat = nullptr;
        this->mpDiffCond = nullptr;
        this->mpEModulus = nullptr;
        this->mpPoisson = nullptr;
        this->mpThermoExpansion = nullptr;
        this->mpYieldStrength = nullptr;
    }

    /**
     * @brief Clone constructor
     */
    Phase::Pointer Phase::Clone() const
    {
        return Kratos::make_shared<Phase>(*this);
    }

    /**
     * @brief create new phase
     */

    Phase::Pointer Phase::Create(Kratos::Parameters NewParameters) const
    {
        if (NewParameters.Has("type"))
        {
            const std::string &type_name = NewParameters["type"].GetString();
            PhaseType type = GetPhaseTypeFromString(type_name); // 字符串转枚举
            return Kratos::make_shared<Phase>(type);
        }
        return Kratos::make_shared<Phase>(PhaseType::UNDEFINED);
    }

    /**
     * Get internal members
     */
    const Phase::Status Phase::Stat()
    {
        /* return the flag of data status */
        return this->mstat;
    }

    const Variable<double>* Phase::GetDensityPtr()
    {
        return (this->mstat == Status::Set ? this->mpDensity : nullptr);
    }

    const Variable<double>* Phase::GetConductivityPtr()
    {
        return (this->mstat == Status::Set ? this->mpConductivity : nullptr);
    }

    const Variable<double>* Phase::GetSpecificHeatPtr()
    {
        return (this->mstat == Status::Set ? this->mpSpecificHeat : nullptr);
    }

    const Variable<double>* Phase::GetDiffCondPtr()
    {
        return (this->mstat == Status::Set ? this->mpDiffCond : nullptr);
    }

    const Variable<double>* Phase::GetEModulusPtr()
    {
        return (this->mstat == Status::Set ? this->mpEModulus : nullptr);
    }

    const Variable<double>* Phase::GetPoissonPtr()
    {
        return (this->mstat == Status::Set ? this->mpPoisson : nullptr);
    }

    const Variable<double>* Phase::GetThermoExpansionPtr()
    {
        return (this->mstat == Status::Set ? this->mpThermoExpansion : nullptr);
    }

    const Variable<double>* Phase::GetYieldStrengthPtr()
    {
        return (this->mstat == Status::Set ? this->mpYieldStrength : nullptr);
    }

    const Phase::TCreepType::Pointer Phase::GetCreepModel()
    {
        return (this->mPhaseType == PhaseType::UNDEFINED ) ? nullptr :this->mpCreep;
    }

    int Phase::Initialize(const Properties& rMaterialProperties, PhaseType rType) 
    {
        int iret = 0x0;
        bool chk;

        this->mPhaseType = rType;

        KRATOS_TRY;

        if(this->mstat == Status::UnSet)
        {
            switch (mPhaseType)
            {
            case PhaseType::AUSTENITE:
                this->mpDensity = &DENSITY_AUSTENITE;
                this->mpConductivity = &CONDUCTIVITY_AUSTENITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_AUSTENITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_AUSTENITE;
                this->mpEModulus = &ELASTICITY_MODULUS_AUSTENITE;
                this->mpPoisson = &POISSON_RATIO_AUSTENITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_AUSTENITE;
                this->mpYieldStrength = &YIELD_STRENGTH_AUSTENITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_AUSTENITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::FERRITE:
                this->mpDensity = &DENSITY_FERRITE;
                this->mpConductivity = &CONDUCTIVITY_FERRITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_FERRITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_FERRITE;
                this->mpEModulus = &ELASTICITY_MODULUS_FERRITE;
                this->mpPoisson = &POISSON_RATIO_FERRITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_FERRITE;
                this->mpYieldStrength = &YIELD_STRENGTH_FERRITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_FERRITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::CEMENTITE:
                this->mpDensity = &DENSITY_CEMENTITE;
                this->mpConductivity = &CONDUCTIVITY_CEMENTITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_CEMENTITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_CEMENTITE;
                this->mpEModulus = &ELASTICITY_MODULUS_CEMENTITE;
                this->mpPoisson = &POISSON_RATIO_CEMENTITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_CEMENTITE;
                this->mpYieldStrength = &YIELD_STRENGTH_CEMENTITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_CEMENTITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::PEARLITE:
                this->mpDensity = &DENSITY_PEARLITE;
                this->mpConductivity = &CONDUCTIVITY_PEARLITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_PEARLITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_PEARLITE;
                this->mpEModulus = &ELASTICITY_MODULUS_PEARLITE;
                this->mpPoisson = &POISSON_RATIO_PEARLITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_PEARLITE;
                this->mpYieldStrength = &YIELD_STRENGTH_PEARLITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_PEARLITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::UPPER_BAINITE:
                this->mpDensity = &DENSITY_UPPER_BAINITE;
                this->mpConductivity = &CONDUCTIVITY_UPPER_BAINITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_UPPER_BAINITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_UPPER_BAINITE;
                this->mpEModulus = &ELASTICITY_MODULUS_UPPER_BAINITE;
                this->mpPoisson = &POISSON_RATIO_UPPER_BAINITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_UPPER_BAINITE;
                this->mpYieldStrength = &YIELD_STRENGTH_UPPER_BAINITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_UPPER_BAINITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::LOWER_BAINITE:
                this->mpDensity = &DENSITY_LOWER_BAINITE;
                this->mpConductivity = &CONDUCTIVITY_LOWER_BAINITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_LOWER_BAINITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_LOWER_BAINITE;
                this->mpEModulus = &ELASTICITY_MODULUS_LOWER_BAINITE;
                this->mpPoisson = &POISSON_RATIO_LOWER_BAINITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_LOWER_BAINITE;
                this->mpYieldStrength = &YIELD_STRENGTH_LOWER_BAINITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_LOWER_BAINITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::MARTENSITE:
                this->mpDensity = &DENSITY_MARTENSITE;
                this->mpConductivity = &CONDUCTIVITY_MARTENSITE;
                this->mpSpecificHeat = &SPECIFIC_HEAT_MARTENSITE;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY_MARTENSITE;
                this->mpEModulus = &ELASTICITY_MODULUS_MARTENSITE;
                this->mpPoisson = &POISSON_RATIO_MARTENSITE;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT_MARTENSITE;
                this->mpYieldStrength = &YIELD_STRENGTH_MARTENSITE;

                if(this->mpCreep == nullptr)
                {
                    this->mpCreep = rMaterialProperties.GetValue(CREEP_MARTENSITE);
                }
                if(this->mpCreep != nullptr)
                {
                    this->mpCreep->Initialize(rMaterialProperties, this->mPhaseType);
                }

                break;
            case PhaseType::NO_PHASE:
                this->mpDensity = &DENSITY;
                this->mpConductivity = &CONDUCTIVITY;
                this->mpSpecificHeat = &SPECIFIC_HEAT;
                this->mpDiffCond = &DIFFUSE_CONDUCTIVITY;
                this->mpEModulus = &ELASTICITY_MODULUS;
                this->mpPoisson = &POISSON_RATIO;
                this->mpThermoExpansion = &THERMAL_EXPANSION_COEFFICIENT;
                this->mpYieldStrength = &YIELD_STRENGTH;
                this->mpCreep = nullptr;
                break;
            default:
                KRATOS_ERROR << "Phase:: Phase type is invalid!" << std::endl;
                break;
            }
            chk = this->mpDensity != nullptr && this->mpConductivity != nullptr && this->mpSpecificHeat != nullptr &&
                      this->mpDiffCond != nullptr && this->mpEModulus != nullptr && this->mpPoisson != nullptr &&
                      this->mpThermoExpansion != nullptr && this->mpYieldStrength != nullptr;
            if(chk)
            {
                this->mstat = Status::Set;
            }
            else
            {
                this->mstat = Status::UnSet;
                KRATOS_ERROR << "Phase:: Initialize failed, some variable is not assigned!" << std::endl;
            }
            iret = 0x1;
        }

        KRATOS_CATCH("");

        return iret;
    }

    double Phase::GetDensity(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;

        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpDensity, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        
        /* return the data */
        return dret;
    }

    double Phase::GetConductivity(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpConductivity, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetSpecificHeat(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpSpecificHeat, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetDensityMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpDensity, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpDensity, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetConductivityMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpConductivity, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpConductivity, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetSpecificHeatMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpSpecificHeat, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpSpecificHeat, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetDC(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpDiffCond, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetDCMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpDiffCond, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpDiffCond, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetElasticityModulus(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpEModulus, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetPoisson(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpPoisson, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetExpansion(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpThermoExpansion, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetElasticityModulusMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpEModulus, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpEModulus, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetPoissonMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpPoisson, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpPoisson, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetExpansionMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double x1, x2, dret = 0;
        if(this->mstat == Status::Set)
        {
            x1 = rMaterialProperties.GetValue(*mpThermoExpansion, rElementGeometry[NodeId], 0, rCurrentProcessInfo);
            x2 = rMaterialProperties.GetValue(*mpThermoExpansion, rElementGeometry[NodeId], 1, rCurrentProcessInfo);
            dret = 0.5 * (x1 + x2);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetYieldStrength(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            dret = rMaterialProperties.GetValue(*mpYieldStrength, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetHardenRate(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mstat == Status::Set)
        {
            auto accessor_yield = rMaterialProperties.GetAccessor(*mpYieldStrength);
            dret = accessor_yield.GetDerivative(EQUIVALENT_STRAIN_PLASTICITY, *mpYieldStrength, rMaterialProperties, rElementGeometry[NodeId], rCurrentProcessInfo);
        }
        /* return the data */
        return dret;
    }

    int Phase::GetYSandHR(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                             const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo,
                             SizeType NodeId, double *ys, double *hr)
    {
        int iret = 0x0;
        double t;
        t = rCurrentProcessInfo[TIME];

        if(ys != nullptr) *ys = DBL_MAX;
        if(hr != nullptr) *hr = 0;

        if(NodeId > 0 && t > 0 && this->mstat == Status::Set)
        {
            *ys = GetYieldStrength(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
            *hr = GetHardenRate(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
            iret = 0x1;
        }
        /* return the data */
        return iret;
    }

    /**
     * 返回creep model的参数列表和参数个数
     * 参数：1.factor_a 2. factor_stress 3. factor_creep 4. factor_temperature 5. factor_time
     * 参数个数 5
     */
    int Phase::GetCreepParameters(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId, SizeType *pnm, double *pdata)
    {
        int iret = 0x0;
        if(this->mpCreep != nullptr)
        {
            iret = this->mpCreep->GetCreepParameters(rMaterialProperties, rElementGeometry,
                                                     rShapeFunctionsValues, rCurrentProcessInfo,
                                                     NodeId, pnm, pdata);
        }
        /* return the data */
        return iret;
    }

    double Phase::GetEqCreepRate(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mpCreep != nullptr)
        {
            dret = this->mpCreep->GetEqCreepRate(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
        }
        /* return the data */
        return dret;
    }

    double Phase::GetEqCreepRateMean(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                            const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                            SizeType NodeId)
    {
        double dret = 0;
        if(this->mpCreep != nullptr)
        {
            dret = this->mpCreep->GetEqCreepRateMean(rMaterialProperties, rElementGeometry,
                                                    rShapeFunctionsValues, rCurrentProcessInfo, NodeId);
        }
        /* return the data */
        return dret;
    }

}