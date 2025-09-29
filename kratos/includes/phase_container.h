#if !defined(KRATOS_PHASE_MODEL)
#define  KRATOS_PHASE_MODEL

#pragma once

// Kratos includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "containers/flags.h"
#include "includes/process_info.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/properties.h"


namespace Kratos {

class PhaseContainer
{
public:
    
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node > GeometryType;


    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PhaseContainer);

    /// Default constructor
    PhaseContainer() = default;


    /// Update internal state from given temperature and time
    void Update(const Properties& rMaterialProperties,
                const GeometryType& rElementGeometry,
                const ProcessInfo& rCurrentProcessInfo);

    void Reset();

    double GetPhaseFraction(const std::string& phase) const;
    
    void UpdateTransformationStrain(const std::string& phase, double increment);
    void UpdateTripStrain(const std::string& phase, double increment);
    void UpdateCreepStrain(const std::string& phase, double increment);


    /// Serialization
    void save(Serializer& rSerializer) const
    {
        rSerializer.save("PhaseFraction", mPhaseFraction);
        rSerializer.save("TransformationStrain", mTransformationStrain);
        rSerializer.save("TripStrain", mTripStrain);
        rSerializer.save("CreepStrain", mCreepStrain);
        rSerializer.save("LatentHeat", mLatentHeat);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("PhaseFraction", mPhaseFraction);
        rSerializer.load("TransformationStrain", mTransformationStrain);
        rSerializer.load("TripStrain", mTripStrain);
        rSerializer.load("CreepStrain", mCreepStrain);
        rSerializer.load("LatentHeat", mLatentHeat);
    }

    std::map<std::string, double> mPhaseFraction;
    std::map<std::string, Matrix> mTransformationStrain;
    std::map<std::string, Matrix> mTripStrain;
    std::map<std::string, Matrix> mCreepStrain;
    std::map<std::string, double> mLatentHeat;

private:
    
};

} // namespace Kratos
#endif