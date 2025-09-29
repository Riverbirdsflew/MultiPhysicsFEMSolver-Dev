//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  Main authors:    whf
//                   



#if !defined(KRATOS_PHASE)
#define  KRATOS_PHASE

/* System includes */

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
#include "creep/maxwell_linear.h"
#include "creep/maxwell_power.h"
#include "phase_transformation_type.h"

namespace Kratos
{

/**
 * Base class of phase object.
 */
/**
 * @brief Base class of phase object.
 * 1. austenite
 * 2. ferrite
 * 3. cementite
 * 4. pearlite
 * 5. upper bainite
 * 6. lower bainite
 * 7. martensite
 */
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) Phase
{

public:

    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node> GeometryType;
    //typedef intrusive_ptr<Phase> Pointer;
    typedef CreepModel TCreepType;

    /* define the status id */
	enum class Status : int {
		UnSet = 0,
		Set = 1,
	};

    /**
     * Counted pointer of phase
     */
    KRATOS_CLASS_POINTER_DEFINITION(Phase);

    /**
     * Constructor
     */
    Phase();

    /**
     * Constructor with model information.
     */
    Phase(PhaseType type, TCreepType::Pointer pCreepModel = nullptr);

    /**
     * Copy Constructor.
     */
    Phase(const Phase& other);
    
    /**
     * Destructor
     */
    virtual ~Phase() = default;

    virtual Phase::Pointer Clone() const;

    virtual Pointer Create(Kratos::Parameters NewParameters) const;

    /**
     * Get internal members
     */
    virtual const Status Stat();
    const Variable<double>* GetDensityPtr();
    const Variable<double>* GetConductivityPtr();
    const Variable<double>* GetSpecificHeatPtr();
    const Variable<double>* GetDiffCondPtr();
    const Variable<double>* GetEModulusPtr();
    const Variable<double>* GetPoissonPtr();
    const Variable<double>* GetThermoExpansionPtr();
    const Variable<double>* GetYieldStrengthPtr();
    virtual const TCreepType::Pointer GetCreepModel();
    
    PhaseType GetPhaseType() const
    {
        return this->mPhaseType;
    }
    // 返回相类型的字符串
    std::string GetPhaseTypeName() const
    {
        switch (mPhaseType)
        {
            case PhaseType::AUSTENITE:      return "AUSTENITE";
            case PhaseType::FERRITE:        return "FERRITE";
            case PhaseType::CEMENTITE:      return "CEMENTITE";
            case PhaseType::PEARLITE:       return "PEARLITE";
            case PhaseType::UPPER_BAINITE:  return "UPPER_BAINITE";
            case PhaseType::LOWER_BAINITE:  return "LOWER_BAINITE";
            case PhaseType::MARTENSITE:     return "MARTENSITE";
            default:             return "Unknown Phase";
        }
    }


    virtual int Initialize(const Properties& rMaterialProperties, PhaseType rType);

    /* return the density value */
	virtual double GetDensity(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the conductivity value */
	virtual double GetConductivity(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the specificheat value */
	virtual double GetSpecificHeat(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

	/* return the integral mean value of density */
	virtual double GetDensityMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the integral mean value of conductivity */
	virtual double GetConductivityMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the integral mean value of specificheat */
	virtual double GetSpecificHeatMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

	/* return the diffusion coefficient value */
	virtual double GetDC(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the integral mean value of the diffusion coefficient */
	virtual double GetDCMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
   
    /* return the elasticity modulus value */
    virtual double GetYoungModulus(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the poisson value */
	virtual double GetPoisson(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the expansion value */
	virtual double GetExpansion(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

	/* return the integral mean value of elasticity modulus */
	virtual double GetYoungModulusMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the integral mean value of poisson */
	virtual double GetPoissonMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the integral mean value of expansion */
	virtual double GetExpansionMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);

	/* return the yield strength */
	virtual double GetYieldStrength(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the hardening rate */
	virtual double GetHardenRate(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the yield strength and hardening rate */
	virtual int GetYSandHR(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId, double *ys, double *hr);

	/* return the creep cofficient data */
	virtual int GetCreepParameters(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId, SizeType *pnm, double *pdata);
	/* return the equivalent creep strain rate */
	virtual double GetEqCreepRate(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);
	/* return the integral mean value of equivalent creep strain rate */
	virtual double GetEqCreepRateMean(const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       SizeType NodeId);


protected:
    
    Status mstat;
    PhaseType mPhaseType;
    CreepModel::Pointer mpCreep;

    Variable<double>* mpDensity;
    Variable<double>* mpConductivity;
    Variable<double>* mpSpecificHeat;
    
    Variable<double>* mpDiffCond; // diffusion coefficient

    Variable<double>* mpEModulus;
    Variable<double>* mpPoisson;
    Variable<double>* mpThermoExpansion; // thermal expansion coefficient
    Variable<double>* mpYieldStrength;

};

}
#endif // KRATOS_PHASE defined