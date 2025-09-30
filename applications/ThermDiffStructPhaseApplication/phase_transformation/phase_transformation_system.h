//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    WHF
//

#pragma once

#if !defined(KRATOS_MATERIAL_PHASE_TRANSFORMATION_SYSTEM)
#define KRATOS_MATERIAL_PHASE_TRANSFORMATION_SYSTEM

/* System includes */
#include <vector>
#include <map>
#include <memory>
#include <array>
#include <string>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "containers/data_value_container.h"
#include "containers/flags.h"
#include "phase_transformation_law.h"
#include "phase.h"
#include "phase_transformation_type.h"
#include "creep/creep_model.h"
#include "kinetics/kinetics_model.h"
#include "trip/trip_model.h"
#include "partition/partition_model.h"



namespace Kratos
{

/**
 * @brief 材料相变过程总控制系统
 * 
 * 这个类是相变过程的核心控制器，负责：
 * 1. 协调相变动力学模型、配分模型和相变塑性模型
 * 2. 计算相分数变化及其对热场、浓度场和结构场的影响
 * 3. 管理七相钢材组织演化（奥氏体、铁素体、渗碳体、珠光体、马氏体、上/下贝氏体）
 * 4. 统一处理四种相变过程（奥氏体转变、共析转变、贝氏体转变、马氏体转变）
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) PhaseTransformationSystem
{
public:
    /// Type definitions
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node> GeometryType;
    typedef Kratos::intrusive_ptr<PhaseTransformationSystem> Pointer;

    static constexpr int NumPhases = static_cast<int> (PhaseType::NUM_PHASES);
    static constexpr int NumTransformations = static_cast<int> (PhaseTransformationType::NUM_TRANSFORMATIONS);

    /* define the status id */
	enum class Status : int {
		UnSet = 0,
		Set = 1,
	};

    /// calculation flags
    struct CalculationFlags {
        bool is_thermo_calculated = true;                           // 热场计算
        bool is_diffusion_calculated = true;                        // 浓度场计算
        bool is_structure_calculated = true;                        // 结构场计算
        bool is_thermal_strain_calculated = true;                   // 热应变计算
        bool is_phase_transformation_calculated = true;             // 相变计算
        bool is_multiphase_defined = true;                          // 多相定义
        bool is_phase_latent_calculated = true;                     // 相变潜热计算
        bool is_phase_transformation_strain_calculated = true;      // 相变应变计算
        bool is_trip_calculated = true;                             // 相变塑性计算
        bool is_creep_calculated = true;                            // 蠕变计算
    };

    /**
     * Constructor
     */
    PhaseTransformationSystem();
    
    /**
     * Destructor
     */
    virtual ~PhaseTransformationSystem() = default;
    
    /**
     * @brief 返回内部成员
     */
    virtual const Status Stat();
    virtual const Phase::Pointer GetPhaseObjectPtr(unsigned int id);
    virtual const PhaseTransformationLaw::Pointer GetPhaseTransObjectPtr(unsigned int id);

    /**
     * @brief 初始化相变系统
     * 设置各种子模型（动力学、配分、TRIP、蠕变模型）
     */
    int Initialize(const Properties& rMaterialProperties,
                   PhaseTransformationType rTransType,
                   PhaseType rPhaseType);
    
    /* microstructure part */
    /* return 7 phase mass fractions */
    int GetMassFractions(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                              SizeType NodeId, std::array<double, NumPhases>& MassFrac);
    int GetLastMassFractions(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                              SizeType NodeId, std::array<double, NumPhases>& LastMassFrac);

	/* return the density value */
	virtual double GetDensity(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the conductivity value */
	virtual double GetConductivity(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the specificheat value */
	virtual double GetSpecificHeat(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);

	/* return the integral mean value of density */
	virtual double GetDensityMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the integral mean value of conductivity */
	virtual double GetConductivityMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the integral mean value of specificheat */
	virtual double GetSpecificHeatMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);

	/* return the diffusion coefficient value */
	virtual double GetDC(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the integral mean value of the diffusion coefficient */
	virtual double GetDCMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);

	/* return the elasticity modulus value */
	virtual double GetElasticityModulus(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the poisson value */
	virtual double GetPoisson(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the expansion value */
	virtual double GetExpansion(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);

	/* return the integral mean value of elasticity modulus */
	virtual double GetElasticityModulusMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the integral mean value of poisson */
	virtual double GetPoissonMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the integral mean value of expansion */
	virtual double GetExpansionMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);

	/* return the yield strength */
	virtual double GetYieldStrength(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the hardening rate */
	virtual double GetHardenRate(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the yield strength and hardening rate */
	virtual int GetYSandHR(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *ys, double *hr);

	/* no-mechanical strain part */
	/* return the integral mean value of thermal strain component increment */
	virtual double GetThermalStrainInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);
	/* return the phase transformation strain component increment */
	virtual double GetPhaseTransStrainInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *mfinc);
	/* return the creep cofficient data */
	virtual int GetCreepParameters(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, unsigned int *pnm, double *pdata);
	/* return the creep strain increment */
	virtual int GetCreepStrainInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double se, double *s, double *cpinc);
	/* return the integral mean value of creep strain increment */
	virtual int GetCreepStrainIncMean(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double se, double *s, double *cpinc);
	/* return the trip increment */
	virtual int GetTripInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *mf, double *mfinc, double *s, double *tpinc);

	/* return the phase transformation latent */
	virtual double GetLatent(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *mfinc);

	/* phase transformation part */
	/* return the austenite mass fraction increment */
	virtual int GetAusteniteMassFrcInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc);
	/* return the eutectoid phase mass fraction increment */
	virtual int GetEutectoidPhaseMassFrcInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc);
	/* return the bainite mass fraction increment */
	virtual int GetBainiteMassFrcInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc);
	/* return the martensite mass fraction increment */
	virtual int GetMartensiteMassFrcInc(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc);
	/* update the incubation, carbon content and mass fraction of all microstructures */
	virtual int UpdateMassFraction(const Properties& rMaterialProperties, GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId, double *cfinc, double *mfinc);

	/* update the carbon content of all microstructures  according to the mean carbon content */
	virtual int UpdateCarbonContent(const Properties& rMaterialProperties, GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues, const ProcessInfo& rCurrentProcessInfo, 
                                 SizeType NodeId);


protected:

    /**
     *mCalculationFlags, 控制参数
     *mpPhases, pointers to Phase objects for each of 7 phase
     *mpPhaseTrans, pointers to PhaseTransformationLaw objects for each of 4 transformation
     */ 
    Status mstat;
    CalculationFlags mCalculationFlags;
    Phase::Pointer mpPhases[7];
    PhaseTransformationLaw::Pointer mpPhaseTrans[4];

};

} // namespace Kratos

#endif // KRATOS_MATERIAL_PHASE_TRANSFORMATION_SYSTEM