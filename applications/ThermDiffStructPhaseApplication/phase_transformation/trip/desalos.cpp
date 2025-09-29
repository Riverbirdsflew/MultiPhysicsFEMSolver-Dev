//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    WHF
//

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "desalos.h"
#include "../therm_diff_struct_phase_application_variables.h"

namespace Kratos
{

    using SizeType = Desalos::SizeType;
    //************************************************************************************
    //************************************************************************************

    Desalos::Desalos()
    {
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mstat = Status::UnSet;
        this->mPara = nullptr;
    }

    /**
     * Copy Constructor.
     */
    Desalos::Desalos(const Desalos &rOther) 
    {
        this->mType = rOther.mType;
        this->mstat = rOther.mstat;
        this->mPara = rOther.mPara;
    }

    /**
     * Destructor.
     */
    Desalos::~Desalos() 
    {
        this->mstat = Status::UnSet;
        this->mType = PhaseTransformationType::UNDEFINED_TRANSFORMATION;
        this->mPara = nullptr;
    }

    std::string Desalos::GetTransTypeName() const 
    {
        switch (mType) {
            case PhaseTransformationType::UNDEFINED_TRANSFORMATION:
                return "UNDEFINED_TRANSFORMATION";
            case PhaseTransformationType::AUSTENIZATION:
                return "AUSTENIZATION";
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                return "EUTECTOID_DECOMPOSITION";
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                return "BAINITIC_TRANSFORMATION";
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                return "MARTENSITIC_TRANSFORMATION";
            default:
                return "UNKNOWN_TRANSFORMATION";
        }
    }

    //************************************************************************************
    //************************************************************************************

    TripModel::Pointer Desalos::Clone() const
    {
        return Kratos::make_shared<Desalos>(*this);
    }

    //************************************************************************************
    //************************************************************************************

    int Desalos::Initialize(const Properties& rMaterialProperties, PhaseTransformationType rType)
    {
        bool chk;
        int iret = 0x0;
        
        this->mType = rType;

        switch(mType)
        {
            case PhaseTransformationType::AUSTENIZATION:
                this->mPara = &TRIP_PARA_AUSTENIZATION_1;
            case PhaseTransformationType::EUTECTOID_DECOMPOSITION:
                this->mPara = &TRIP_PARA_EUTECTOID_DECOMPOSITION_1;
                break;
            case PhaseTransformationType::BAINITIC_TRANSFORMATION:
                this->mPara = &TRIP_PARA_BAINITIC_TRANSFORMATION_1;
                break;
            case PhaseTransformationType::MARTENSITIC_TRANSFORMATION:
                this->mPara = &TRIP_PARA_MARTENSITIC_TRANSFORMATION_1;
            default:
                KRATOS_ERROR <<  "The Phase Transformation Type is not defined."<< std::endl;
        }

        chk = this->mPara != nullptr;
        if (!chk)
            KRATOS_ERROR <<  "The parameters for Desalos model are not assigned"<< std::endl;
        else
            this->mstat = Status::Set;

        iret = 0x1;

        return iret;

    }

    
    /**
     * 计算相变诱发的TRIP增量
     * 共析分解TRIP参数： 1.factor_k
     */
    // eutectoid trip
    double Desalos::GetEqTripInc(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                    const Vector &rShapeFunctionsValues, const ProcessInfo &rCurrentProcessInfo,
                                    SizeType NodeId, double rMassFrac, double rMassFracInc)
    {

        double k, dt, t;
        double eq_trip_inc = 0;

        dt = rCurrentProcessInfo[DELTA_TIME];
        t = rCurrentProcessInfo[TIME];

        if (t >= 0 && dt >= 0 && rMassFrac >= 0 && rMassFracInc > 0)
        {
            k = rMaterialProperties.GetValue(*mPara, rElementGeometry[NodeId], rCurrentProcessInfo);
            eq_trip_inc = 3 * k * (1 - rMassFrac) * rMassFracInc;
        }
        return eq_trip_inc;
    }
}
