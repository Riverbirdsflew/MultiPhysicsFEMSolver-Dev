//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pablo Becker
//

#if !defined(KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED )
#define  KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/define.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Convection diffusion settings. This class contains information to be used by the convection diffusion solver, all the variables that will be needed by the solver.
/** All the variables needed by any convection diffusion problem are included here. However, many problems do not require such a large amount of variables.
 * For those cases, variables that are not defined will be set to either zero or one (depending on the variable)
 * For this purpose, there are flags to ask whether a variable has been defined or not. For each variable, there are three main functions.
 * SetVariableforXuse: we assign the variable for this use. When doing that we also set the flag that now this variable has been defined.
 * GetVariableforXuse: we return the variable for this use.
 * IsDefinedVariableforXuse: tells whether that variable has been defined or not.
*/
class KRATOS_API(KRATOS_CORE) ConvectionDiffusionSettings
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConvectionDiffusionSettings
    KRATOS_CLASS_POINTER_DEFINITION(ConvectionDiffusionSettings);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConvectionDiffusionSettings() = default;

    ConvectionDiffusionSettings(const ConvectionDiffusionSettings& rOther):
        mpDensityVar(rOther.mpDensityVar),
        mpDensityM1Var(rOther.mpDensityM1Var),//多相密度7
        mpDensityM2Var(rOther.mpDensityM2Var),
        mpDensityM3Var(rOther.mpDensityM3Var),
        mpDensityM4Var(rOther.mpDensityM4Var),
        mpDensityM5Var(rOther.mpDensityM5Var),
        mpDensityM6Var(rOther.mpDensityM6Var),
        mpDensityM7Var(rOther.mpDensityM7Var),
        mpMassFractionM1Var(rOther.mpMassFractionM1Var),//多相质量分数7
        mpMassFractionM2Var(rOther.mpMassFractionM2Var),
        mpMassFractionM3Var(rOther.mpMassFractionM3Var),
        mpMassFractionM4Var(rOther.mpMassFractionM4Var),
        mpMassFractionM5Var(rOther.mpMassFractionM5Var),
        mpMassFractionM6Var(rOther.mpMassFractionM6Var),
        mpMassFractionM7Var(rOther.mpMassFractionM7Var),
        mpDiffusionVar(rOther.mpDiffusionVar),
        mpDiffusionM1Var(rOther.mpDiffusionM1Var),//多相导热系数7
        mpDiffusionM2Var(rOther.mpDiffusionM2Var),
        mpDiffusionM3Var(rOther.mpDiffusionM3Var),
        mpDiffusionM4Var(rOther.mpDiffusionM4Var),
        mpDiffusionM5Var(rOther.mpDiffusionM5Var),
        mpDiffusionM6Var(rOther.mpDiffusionM6Var),
        mpDiffusionM7Var(rOther.mpDiffusionM7Var),
        mpUnknownVar(rOther.mpUnknownVar),
        mpVolumeSourceVar(rOther.mpVolumeSourceVar),
        mpAustenizeLatentVar(rOther.mpAustenizeLatentVar),//相变潜热
        mpEutectoidLatentVar(rOther.mpEutectoidLatentVar),
        mpBainiteLatentVar(rOther.mpEutectoidLatentVar),
        mpMartensiteLatentVar(rOther.mpEutectoidLatentVar),
        mpSurfaceSourceVar(rOther.mpSurfaceSourceVar),
        mpProjectionVar(rOther. mpProjectionVar),
        mpConvectionVar(rOther.mpConvectionVar),
        mpGradientVar(rOther.mpGradientVar),
        mpMeshVelocityVar(rOther.mpMeshVelocityVar),
        mpTransferCoefficientVar(rOther.mpTransferCoefficientVar),
        mpVelocityVar(rOther.mpVelocityVar),
        mpSpecificHeatVar(rOther.mpSpecificHeatVar),
        mpSpecificHeatM1Var(rOther.mpSpecificHeatM1Var),//多相比热7
        mpSpecificHeatM2Var(rOther.mpSpecificHeatM2Var),
        mpSpecificHeatM3Var(rOther.mpSpecificHeatM3Var),
        mpSpecificHeatM4Var(rOther.mpSpecificHeatM4Var),
        mpSpecificHeatM5Var(rOther.mpSpecificHeatM5Var),
        mpSpecificHeatM6Var(rOther.mpSpecificHeatM6Var),
        mpSpecificHeatM7Var(rOther.mpSpecificHeatM7Var),
        mpReactionVar(rOther.mpReactionVar),
        mpReactionGradientVar(rOther.mpReactionGradientVar),
        mis_defined_DensityVar(rOther.mis_defined_DensityVar),
        mis_defined_DensityM1Var(rOther.mis_defined_DensityM1Var),//多相密度
        mis_defined_DensityM2Var(rOther.mis_defined_DensityM2Var),
        mis_defined_DensityM3Var(rOther.mis_defined_DensityM3Var),
        mis_defined_DensityM4Var(rOther.mis_defined_DensityM4Var),
        mis_defined_DensityM5Var(rOther.mis_defined_DensityM5Var),
        mis_defined_DensityM6Var(rOther.mis_defined_DensityM6Var),
        mis_defined_DensityM7Var(rOther.mis_defined_DensityM7Var),
		mis_defined_DiffusionVar(rOther.mis_defined_DiffusionVar),
        mis_defined_DiffusionM1Var(rOther.mis_defined_DiffusionM1Var),//多相热导
        mis_defined_DiffusionM2Var(rOther.mis_defined_DiffusionM2Var),
        mis_defined_DiffusionM3Var(rOther.mis_defined_DiffusionM3Var),
        mis_defined_DiffusionM4Var(rOther.mis_defined_DiffusionM4Var),
        mis_defined_DiffusionM5Var(rOther.mis_defined_DiffusionM5Var),
        mis_defined_DiffusionM6Var(rOther.mis_defined_DiffusionM6Var),
        mis_defined_DiffusionM7Var(rOther.mis_defined_DiffusionM7Var),
		mis_defined_UnknownVar(rOther.mis_defined_UnknownVar),
		mis_defined_VolumeSourceVar(rOther.mis_defined_VolumeSourceVar),
		mis_defined_AustenizeLatentVar(rOther.mis_defined_AustenizeLatentVar),
		mis_defined_EutectoidLatentVar(rOther.mis_defined_EutectoidLatentVar),
		mis_defined_BainiteLatentVar(rOther.mis_defined_BainiteLatentVar),
		mis_defined_MartensiteLatentVar(rOther.mis_defined_MartensiteLatentVar),
		mis_defined_SurfaceSourceVar(rOther.mis_defined_SurfaceSourceVar),
		mis_defined_ProjectionVar(rOther.mis_defined_ProjectionVar),
		mis_defined_ConvectionVar(rOther.mis_defined_ConvectionVar),
        mis_defined_GradientVar(rOther.mis_defined_GradientVar),
		mis_defined_MeshVelocityVar(rOther.mis_defined_MeshVelocityVar),
		mis_defined_TransferCoefficientVar(rOther.mis_defined_TransferCoefficientVar),
		mis_defined_VelocityVar(rOther.mis_defined_VelocityVar),
		mis_defined_SpecificHeatVar(rOther.mis_defined_SpecificHeatVar),
        mis_defined_SpecificHeatM1Var(rOther.mis_defined_SpecificHeatM1Var),//多相比热
        mis_defined_SpecificHeatM2Var(rOther.mis_defined_SpecificHeatM2Var),
        mis_defined_SpecificHeatM3Var(rOther.mis_defined_SpecificHeatM3Var),
        mis_defined_SpecificHeatM4Var(rOther.mis_defined_SpecificHeatM4Var),
        mis_defined_SpecificHeatM5Var(rOther.mis_defined_SpecificHeatM5Var),
        mis_defined_SpecificHeatM6Var(rOther.mis_defined_SpecificHeatM6Var),
        mis_defined_SpecificHeatM7Var(rOther.mis_defined_SpecificHeatM7Var),
        mis_defined_ReactionVar(rOther.mis_defined_ReactionVar),
        mIsDefinedReactionGradientVar(rOther.mIsDefinedReactionGradientVar)
    {
    }

    /// Destructor.
    virtual ~ConvectionDiffusionSettings() {};

    ///@}
    ///@name Operators
    ///@{
    void SetDensityVariable(const Variable<double>& rvar)
    {
        mpDensityVar = &rvar;
        mis_defined_DensityVar=true;
    }
    const Variable<double>& GetDensityVariable() const
    {
        return *mpDensityVar;
    }
    bool IsDefinedDensityVariable() const
    {
		return mpDensityVar != nullptr;
	}
    //多相材料密度
    //m1
    void SetDensityM1Variable(const Variable<double>& rvar)
    {
        mpDensityM1Var = &rvar;
        mis_defined_DensityM1Var=true;
    }
    const Variable<double>& GetDensityM1Variable() const
    {
        return *mpDensityM1Var;
    }
    bool IsDefinedDensityM1Variable() const
    {
		return mpDensityM1Var != nullptr;
	}
    
    void SetMassFractionM1Variable(const Variable<double>& rvar)
    {
        mpMassFractionM1Var = &rvar;
        mis_defined_MassFractionM1Var=true;
    }
    const Variable<double>& GetMassFractionM1Variable() const
    {
        return *mpMassFractionM1Var;
    }
    bool IsDefinedMassFractionM1Variable() const
    {
		return mpMassFractionM1Var != nullptr;
	}
    //m2
    void SetDensityM2Variable(const Variable<double>& rvar)
    {
        mpDensityM2Var = &rvar;
        mis_defined_DensityM2Var=true;
    }
    const Variable<double>& GetDensityM2Variable() const
    {
        return *mpDensityM2Var;
    }
    bool IsDefinedDensityM2Variable() const
    {
		return mpDensityM2Var != nullptr;
	}
    
    void SetMassFractionM2Variable(const Variable<double>& rvar)
    {
        mpMassFractionM2Var = &rvar;
        mis_defined_MassFractionM2Var=true;
    }
    const Variable<double>& GetMassFractionM2Variable() const
    {
        return *mpMassFractionM2Var;
    }
    bool IsDefinedMassFractionM2Variable() const
    {
		return mpMassFractionM2Var != nullptr;
	}
    //m3
    void SetDensityM3Variable(const Variable<double>& rvar)
    {
        mpDensityM3Var = &rvar;
        mis_defined_DensityM3Var=true;
    }
    const Variable<double>& GetDensityM3Variable() const
    {
        return *mpDensityM3Var;
    }
    bool IsDefinedDensityM3Variable() const
    {
		return mpDensityM3Var != nullptr;
	}
    
    void SetMassFractionM3Variable(const Variable<double>& rvar)
    {
        mpMassFractionM3Var = &rvar;
        mis_defined_MassFractionM3Var=true;
    }
    const Variable<double>& GetMassFractionM3Variable() const
    {
        return *mpMassFractionM3Var;
    }
    bool IsDefinedMassFractionM3Variable() const
    {
		return mpMassFractionM3Var != nullptr;
	}
    //m4
    void SetDensityM4Variable(const Variable<double>& rvar)
    {
        mpDensityM4Var = &rvar;
        mis_defined_DensityM4Var=true;
    }
    const Variable<double>& GetDensityM4Variable() const
    {
        return *mpDensityM4Var;
    }
    bool IsDefinedDensityM4Variable() const
    {
		return mpDensityM4Var != nullptr;
	}
    
    void SetMassFractionM4Variable(const Variable<double>& rvar)
    {
        mpMassFractionM4Var = &rvar;
        mis_defined_MassFractionM4Var=true;
    }
    const Variable<double>& GetMassFractionM4Variable() const
    {
        return *mpMassFractionM4Var;
    }
    bool IsDefinedMassFractionM4Variable() const
    {
		return mpMassFractionM4Var != nullptr;
	}
    //m5
    void SetDensityM5Variable(const Variable<double>& rvar)
    {
        mpDensityM5Var = &rvar;
        mis_defined_DensityM5Var=true;
    }
    const Variable<double>& GetDensityM5Variable() const
    {
        return *mpDensityM5Var;
    }
    bool IsDefinedDensityM5Variable() const
    {
		return mpDensityM5Var != nullptr;
	}
    
    void SetMassFractionM5Variable(const Variable<double>& rvar)
    {
        mpMassFractionM5Var = &rvar;
        mis_defined_MassFractionM5Var=true;
    }
    const Variable<double>& GetMassFractionM5Variable() const
    {
        return *mpMassFractionM5Var;
    }
    bool IsDefinedMassFractionM5Variable() const
    {
		return mpMassFractionM5Var != nullptr;
	}
    //m6
    void SetDensityM6Variable(const Variable<double>& rvar)
    {
        mpDensityM6Var = &rvar;
        mis_defined_DensityM6Var=true;
    }
    const Variable<double>& GetDensityM6Variable() const
    {
        return *mpDensityM6Var;
    }
    bool IsDefinedDensityM6Variable() const
    {
		return mpDensityM6Var != nullptr;
	}
    
    void SetMassFractionM6Variable(const Variable<double>& rvar)
    {
        mpMassFractionM6Var = &rvar;
        mis_defined_MassFractionM6Var=true;
    }
    const Variable<double>& GetMassFractionM6Variable() const
    {
        return *mpMassFractionM6Var;
    }
    bool IsDefinedMassFractionM6Variable() const
    {
		return mpMassFractionM6Var != nullptr;
	}
    //m7
    void SetDensityM7Variable(const Variable<double>& rvar)
    {
        mpDensityM7Var = &rvar;
        mis_defined_DensityM7Var=true;
    }
    const Variable<double>& GetDensityM7Variable() const
    {
        return *mpDensityM7Var;
    }
    bool IsDefinedDensityM7Variable() const
    {
		return mpDensityM7Var != nullptr;
	}
    
    void SetMassFractionM7Variable(const Variable<double>& rvar)
    {
        mpMassFractionM7Var = &rvar;
        mis_defined_MassFractionM7Var=true;
    }
    const Variable<double>& GetMassFractionM7Variable() const
    {
        return *mpMassFractionM7Var;
    }
    bool IsDefinedMassFractionM7Variable() const
    {
		return mpMassFractionM7Var != nullptr;
	}
    
    
    void SetDiffusionVariable(const Variable<double>& rvar)
    {
        mpDiffusionVar = &rvar;
		mis_defined_DiffusionVar=true;
    }
    const Variable<double>& GetDiffusionVariable() const
    {
        return *mpDiffusionVar;
    }
    bool IsDefinedDiffusionVariable() const
    {
		return mpDiffusionVar != nullptr;
	}
    //多相材料热导
    //m1
    void SetDiffusionM1Variable(const Variable<double>& rvar)
    {
        mpDiffusionM1Var = &rvar;
        mis_defined_DiffusionM1Var=true;
    }
    const Variable<double>& GetDiffusionM1Variable() const
    {
        return *mpDiffusionM1Var;
    }
    bool IsDefinedDiffusionM1Variable() const
    {
		return mpDiffusionM1Var != nullptr;
	}
    
    //m2
    void SetDiffusionM2Variable(const Variable<double>& rvar)
    {
        mpDiffusionM2Var = &rvar;
        mis_defined_DiffusionM2Var=true;
    }
    const Variable<double>& GetDiffusionM2Variable() const
    {
        return *mpDiffusionM2Var;
    }
    bool IsDefinedDiffusionM2Variable() const
    {
		return mpDiffusionM2Var != nullptr;
	}
    
    //m3
    void SetDiffusionM3Variable(const Variable<double>& rvar)
    {
        mpDiffusionM3Var = &rvar;
        mis_defined_DiffusionM3Var=true;
    }
    const Variable<double>& GetDiffusionM3Variable() const
    {
        return *mpDiffusionM3Var;
    }
    bool IsDefinedDiffusionM3Variable() const
    {
		return mpDiffusionM3Var != nullptr;
	}
    
    //m4
    void SetDiffusionM4Variable(const Variable<double>& rvar)
    {
        mpDiffusionM4Var = &rvar;
        mis_defined_DiffusionM4Var=true;
    }
    const Variable<double>& GetDiffusionM4Variable() const
    {
        return *mpDiffusionM4Var;
    }
    bool IsDefinedDiffusionM4Variable() const
    {
		return mpDiffusionM4Var != nullptr;
	}
    
    //m5
    void SetDiffusionM5Variable(const Variable<double>& rvar)
    {
        mpDiffusionM5Var = &rvar;
        mis_defined_DiffusionM5Var=true;
    }
    const Variable<double>& GetDiffusionM5Variable() const
    {
        return *mpDiffusionM5Var;
    }
    bool IsDefinedDiffusionM5Variable() const
    {
		return mpDiffusionM5Var != nullptr;
	}
    
    //m6
    void SetDiffusionM6Variable(const Variable<double>& rvar)
    {
        mpDiffusionM6Var = &rvar;
        mis_defined_DiffusionM6Var=true;
    }
    const Variable<double>& GetDiffusionM6Variable() const
    {
        return *mpDiffusionM6Var;
    }
    bool IsDefinedDiffusionM6Variable() const
    {
		return mpDiffusionM6Var != nullptr;
	}
    
    //m7
    void SetDiffusionM7Variable(const Variable<double>& rvar)
    {
        mpDiffusionM7Var = &rvar;
        mis_defined_DiffusionM7Var=true;
    }
    const Variable<double>& GetDiffusionM7Variable() const
    {
        return *mpDiffusionM7Var;
    }
    bool IsDefinedDiffusionM7Variable() const
    {
		return mpDiffusionM7Var != nullptr;
	}
    

    void SetUnknownVariable(const Variable<double>& rvar)
    {
        mpUnknownVar = &rvar;
		mis_defined_UnknownVar=true;
    }
    const Variable<double>& GetUnknownVariable() const
    {
        return *mpUnknownVar;
    }
    bool IsDefinedUnknownVariable() const
    {
		return mpUnknownVar != nullptr;
	}

    void SetVolumeSourceVariable(const Variable<double>& rvar)
    {
        mpVolumeSourceVar = &rvar;
		mis_defined_VolumeSourceVar=true;
    }
    const Variable<double>& GetVolumeSourceVariable() const
    {
        return *mpVolumeSourceVar;
    }
    bool IsDefinedVolumeSourceVariable() const
    {
		return mpVolumeSourceVar != nullptr;
	}
	
    void SetLatentAustenizeVariable(const Variable<double>& rvar)
    {
        mpAustenizeLatentVar = &rvar;
		mis_defined_AustenizeLatentVar=true;
    }
    const Variable<double>& GetLatentAustenizeVariable() const
    {
        return *mpAustenizeLatentVar;
    }
    bool IsDefinedLatentAustenizeVariable() const
    {
		return mpAustenizeLatentVar != nullptr;
	}
    
    void SetLatentEutectoidVariable(const Variable<double>& rvar)
    {
        mpEutectoidLatentVar = &rvar;
		mis_defined_EutectoidLatentVar=true;
    }
    const Variable<double>& GetLatentEutectoidVariable() const
    {
        return *mpEutectoidLatentVar;
    }
    bool IsDefinedLatentEutectoidVariable() const
    {
		return mpEutectoidLatentVar != nullptr;
	}
	
    void SetLatentBainiteVariable(const Variable<double>& rvar)
    {
        mpBainiteLatentVar = &rvar;
		mis_defined_BainiteLatentVar=true;
    }
    const Variable<double>& GetLatentBainiteVariable() const
    {
        return *mpBainiteLatentVar;
    }
    bool IsDefinedLatentBainiteVariable() const
    {
		return mpBainiteLatentVar != nullptr;
	}

    void SetLatentMartensiteVariable(const Variable<double>& rvar)
    {
        mpMartensiteLatentVar = &rvar;
		mis_defined_MartensiteLatentVar=true;
    }
    const Variable<double>& GetLatentMartensiteVariable() const
    {
        return *mpMartensiteLatentVar;
    }
    bool IsDefinedLatentMartensiteVariable() const
    {
		return mpMartensiteLatentVar != nullptr;
	}

    void SetSurfaceSourceVariable(const Variable<double>& rvar)
    {
        mpSurfaceSourceVar = &rvar;
		mis_defined_SurfaceSourceVar=true;
    }
    const Variable<double>& GetSurfaceSourceVariable() const
    {
        return *mpSurfaceSourceVar;
    }
    bool IsDefinedSurfaceSourceVariable() const
    {
		return mpSurfaceSourceVar != nullptr;
	}

    void SetProjectionVariable(const Variable<double>& rvar)
    {
        mpProjectionVar = &rvar;
		mis_defined_ProjectionVar=true;
    }
    const Variable<double>& GetProjectionVariable() const
    {
        return *mpProjectionVar;
    }
    bool IsDefinedProjectionVariable() const
    {
		return mpProjectionVar != nullptr;
	}

    void SetConvectionVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpConvectionVar = &rvar;
		mis_defined_ConvectionVar=true;
    }
    const Variable<array_1d<double,3> >& GetConvectionVariable() const
    {
        return *mpConvectionVar;
    }
    bool IsDefinedConvectionVariable() const
    {
		return mpConvectionVar != nullptr;
	}

    void SetGradientVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpGradientVar = &rvar;
        mis_defined_GradientVar=true;
    }
    const Variable<array_1d<double,3> >& GetGradientVariable() const
    {
        return *mpGradientVar;
    }
    bool IsDefinedGradientVariable() const
    {
        return mpGradientVar != nullptr;
    }

    void SetMeshVelocityVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpMeshVelocityVar = &rvar;
		mis_defined_MeshVelocityVar=true;
    }
    const Variable<array_1d<double,3> >& GetMeshVelocityVariable() const
    {
        return *mpMeshVelocityVar;
    }
    bool IsDefinedMeshVelocityVariable() const
    {
		return mpMeshVelocityVar != nullptr;
	}

    void SetTransferCoefficientVariable(const Variable<double>& rvar)
    {
        mpTransferCoefficientVar = &rvar;
		mis_defined_TransferCoefficientVar=true;
    }
    const Variable<double>& GetTransferCoefficientVariable() const
    {
        return *mpTransferCoefficientVar;
    }
    bool IsDefinedTransferCoefficientVariable() const
    {
		return mpTransferCoefficientVar != nullptr;
	}

    void SetVelocityVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpVelocityVar = &rvar;
		mis_defined_VelocityVar=true;
    }
    const Variable<array_1d<double,3> >& GetVelocityVariable() const
    {
        return *mpVelocityVar;
    }
    bool IsDefinedVelocityVariable() const
    {
		return mpVelocityVar != nullptr;
	}

    void SetSpecificHeatVariable(const Variable<double>& rvar)
    {
        mpSpecificHeatVar = &rvar;
		mis_defined_SpecificHeatVar=true;
    }
    const Variable<double>& GetSpecificHeatVariable() const
    {
        return *mpSpecificHeatVar;
    }
    bool IsDefinedSpecificHeatVariable() const
    {
		return mpSpecificHeatVar != nullptr;
	}
    //多相材料比热
    //m1
    void SetSpecificHeatM1Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM1Var = &rvar;
        mis_defined_SpecificHeatM1Var=true;
    }
    const Variable<double>& GetSpecificHeatM1Variable() const
    {
        return *mpSpecificHeatM1Var;
    }
    bool IsDefinedSpecificHeatM1Variable() const
    {
		return mpSpecificHeatM1Var != nullptr;
	}
    
    //m2
    void SetSpecificHeatM2Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM2Var = &rvar;
        mis_defined_SpecificHeatM2Var=true;
    }
    const Variable<double>& GetSpecificHeatM2Variable() const
    {
        return *mpSpecificHeatM2Var;
    }
    bool IsDefinedSpecificHeatM2Variable() const
    {
		return mpSpecificHeatM2Var != nullptr;
	}
    
    //m3
    void SetSpecificHeatM3Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM3Var = &rvar;
        mis_defined_SpecificHeatM3Var=true;
    }
    const Variable<double>& GetSpecificHeatM3Variable() const
    {
        return *mpSpecificHeatM3Var;
    }
    bool IsDefinedSpecificHeatM3Variable() const
    {
		return mpSpecificHeatM3Var != nullptr;
	}
    
    //m4
    void SetSpecificHeatM4Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM4Var = &rvar;
        mis_defined_SpecificHeatM4Var=true;
    }
    const Variable<double>& GetSpecificHeatM4Variable() const
    {
        return *mpSpecificHeatM4Var;
    }
    bool IsDefinedSpecificHeatM4Variable() const
    {
		return mpSpecificHeatM4Var != nullptr;
	}
    
    //m5
    void SetSpecificHeatM5Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM5Var = &rvar;
        mis_defined_SpecificHeatM5Var=true;
    }
    const Variable<double>& GetSpecificHeatM5Variable() const
    {
        return *mpSpecificHeatM5Var;
    }
    bool IsDefinedSpecificHeatM5Variable() const
    {
		return mpSpecificHeatM5Var != nullptr;
	}
    
    //m6
    void SetSpecificHeatM6Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM6Var = &rvar;
        mis_defined_SpecificHeatM6Var=true;
    }
    const Variable<double>& GetSpecificHeatM6Variable() const
    {
        return *mpSpecificHeatM6Var;
    }
    bool IsDefinedSpecificHeatM6Variable() const
    {
		return mpSpecificHeatM6Var != nullptr;
	}
    
    //m7
    void SetSpecificHeatM7Variable(const Variable<double>& rvar)
    {
        mpSpecificHeatM7Var = &rvar;
        mis_defined_SpecificHeatM7Var=true;
    }
    const Variable<double>& GetSpecificHeatM7Variable() const
    {
        return *mpSpecificHeatM7Var;
    }
    bool IsDefinedSpecificHeatM7Variable() const
    {
		return mpSpecificHeatM7Var != nullptr;
	}

    void SetReactionVariable(const Variable<double>& rvar)
    {
        mpReactionVar = &rvar;
		mis_defined_ReactionVar=true;
    }
    const Variable<double>& GetReactionVariable() const
    {
        return *mpReactionVar;
    }
    bool IsDefinedReactionVariable() const
    {
		return mpReactionVar != nullptr;
	}

    void SetReactionGradientVariable(const Variable<array_1d<double,3>>& rVar)
    {
        mpReactionGradientVar = &rVar;
		mIsDefinedReactionGradientVar=true;
    }
    const Variable<array_1d<double,3>>& GetReactionGradientVariable() const
    {
        return *mpReactionGradientVar;
    }
    bool IsDefinedReactionGradientVariable() const
    {
		return mpReactionGradientVar != nullptr;
	}

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{
    /// Assignment operator.
    ConvectionDiffusionSettings& operator=(ConvectionDiffusionSettings const& rOther)
    {
        mpDensityVar = rOther.mpDensityVar;
        mpDensityM1Var = rOther.mpDensityM1Var;//多相
        mpDensityM2Var = rOther.mpDensityM2Var;
        mpDensityM3Var = rOther.mpDensityM3Var;
        mpDensityM4Var = rOther.mpDensityM4Var;
        mpDensityM5Var = rOther.mpDensityM5Var;
        mpDensityM6Var = rOther.mpDensityM6Var;
        mpDensityM7Var = rOther.mpDensityM7Var;
        mpMassFractionM1Var = rOther.mpMassFractionM1Var;
        mpMassFractionM2Var = rOther.mpMassFractionM2Var;
        mpMassFractionM3Var = rOther.mpMassFractionM3Var;
        mpMassFractionM4Var = rOther.mpMassFractionM4Var;
        mpMassFractionM5Var = rOther.mpMassFractionM5Var;
        mpMassFractionM6Var = rOther.mpMassFractionM6Var;
        mpMassFractionM7Var = rOther.mpMassFractionM7Var;//多相
        mpDiffusionVar = rOther.mpDiffusionVar;
        mpDiffusionM1Var = rOther.mpDiffusionM1Var;//热导
        mpDiffusionM2Var = rOther.mpDiffusionM2Var;
        mpDiffusionM3Var = rOther.mpDiffusionM3Var;
        mpDiffusionM4Var = rOther.mpDiffusionM4Var;
        mpDiffusionM5Var = rOther.mpDiffusionM5Var;
        mpDiffusionM6Var = rOther.mpDiffusionM6Var;
        mpDiffusionM7Var = rOther.mpDiffusionM7Var;
        mpUnknownVar = rOther.mpUnknownVar;
        mpVolumeSourceVar = rOther.mpVolumeSourceVar;
        mpAustenizeLatentVar = rOther.mpAustenizeLatentVar;
        mpEutectoidLatentVar = rOther.mpEutectoidLatentVar;
        mpBainiteLatentVar = rOther.mpEutectoidLatentVar;
        mpMartensiteLatentVar = rOther.mpEutectoidLatentVar;
        mpSurfaceSourceVar = rOther.mpSurfaceSourceVar;
        mpProjectionVar = rOther.mpProjectionVar;
        mpConvectionVar = rOther.mpConvectionVar;
        mpGradientVar = rOther.mpGradientVar;
        mpMeshVelocityVar = rOther.mpMeshVelocityVar;
        mpTransferCoefficientVar = rOther.mpTransferCoefficientVar;
        mpVelocityVar = rOther.mpVelocityVar;
		mpSpecificHeatVar = rOther.mpSpecificHeatVar;
        mpSpecificHeatM1Var = rOther.mpSpecificHeatM1Var;//比热
        mpSpecificHeatM2Var = rOther.mpSpecificHeatM2Var;
        mpSpecificHeatM3Var = rOther.mpSpecificHeatM3Var;
        mpSpecificHeatM4Var = rOther.mpSpecificHeatM4Var;
        mpSpecificHeatM5Var = rOther.mpSpecificHeatM5Var;
        mpSpecificHeatM6Var = rOther.mpSpecificHeatM6Var;
        mpSpecificHeatM7Var = rOther.mpSpecificHeatM7Var;
        mpReactionVar = rOther.mpReactionVar;
        mpReactionGradientVar = rOther.mpReactionGradientVar;
        //now the is_defined
        mis_defined_DensityVar = rOther.mis_defined_DensityVar;
        mis_defined_DensityM1Var = rOther.mis_defined_DensityM1Var;//多相
        mis_defined_DensityM2Var = rOther.mis_defined_DensityM2Var;
        mis_defined_DensityM3Var = rOther.mis_defined_DensityM3Var;
        mis_defined_DensityM4Var = rOther.mis_defined_DensityM4Var;
        mis_defined_DensityM5Var = rOther.mis_defined_DensityM5Var;
        mis_defined_DensityM6Var = rOther.mis_defined_DensityM6Var;
        mis_defined_DensityM7Var = rOther.mis_defined_DensityM7Var;
        mis_defined_MassFractionM1Var = rOther.mis_defined_MassFractionM1Var;
        mis_defined_MassFractionM2Var = rOther.mis_defined_MassFractionM2Var;
        mis_defined_MassFractionM3Var = rOther.mis_defined_MassFractionM3Var;
        mis_defined_MassFractionM4Var = rOther.mis_defined_MassFractionM4Var;
        mis_defined_MassFractionM5Var = rOther.mis_defined_MassFractionM5Var;
        mis_defined_MassFractionM6Var = rOther.mis_defined_MassFractionM6Var;
        mis_defined_MassFractionM7Var = rOther.mis_defined_MassFractionM7Var;//多相
		mis_defined_DiffusionVar = rOther.mis_defined_DiffusionVar;
        mis_defined_DiffusionM1Var = rOther.mis_defined_DiffusionM1Var;//导热
        mis_defined_DiffusionM2Var = rOther.mis_defined_DiffusionM2Var;
        mis_defined_DiffusionM3Var = rOther.mis_defined_DiffusionM3Var;
        mis_defined_DiffusionM4Var = rOther.mis_defined_DiffusionM4Var;
        mis_defined_DiffusionM5Var = rOther.mis_defined_DiffusionM5Var;
        mis_defined_DiffusionM6Var = rOther.mis_defined_DiffusionM6Var;
        mis_defined_DiffusionM7Var = rOther.mis_defined_DiffusionM7Var;
		mis_defined_UnknownVar = rOther.mis_defined_UnknownVar;
		mis_defined_VolumeSourceVar = rOther.mis_defined_VolumeSourceVar;
		mis_defined_AustenizeLatentVar = rOther.mis_defined_AustenizeLatentVar;
                mis_defined_EutectoidLatentVar = rOther.mis_defined_EutectoidLatentVar;
                mis_defined_BainiteLatentVar = rOther.mis_defined_BainiteLatentVar;
                mis_defined_MartensiteLatentVar = rOther.mis_defined_MartensiteLatentVar;
		mis_defined_SurfaceSourceVar = rOther.mis_defined_SurfaceSourceVar;
		mis_defined_ProjectionVar = rOther.mis_defined_ProjectionVar;
        mis_defined_ConvectionVar = rOther.mis_defined_ConvectionVar;
        mis_defined_GradientVar = rOther.mis_defined_GradientVar;
		mis_defined_MeshVelocityVar = rOther.mis_defined_MeshVelocityVar;
		mis_defined_TransferCoefficientVar = rOther.mis_defined_TransferCoefficientVar;
		mis_defined_VelocityVar = rOther.mis_defined_VelocityVar;
		mis_defined_SpecificHeatVar = rOther.mis_defined_SpecificHeatVar;
        mis_defined_SpecificHeatM1Var = rOther.mis_defined_SpecificHeatM1Var;//比热
        mis_defined_SpecificHeatM2Var = rOther.mis_defined_SpecificHeatM2Var;
        mis_defined_SpecificHeatM3Var = rOther.mis_defined_SpecificHeatM3Var;
        mis_defined_SpecificHeatM4Var = rOther.mis_defined_SpecificHeatM4Var;
        mis_defined_SpecificHeatM5Var = rOther.mis_defined_SpecificHeatM5Var;
        mis_defined_SpecificHeatM6Var = rOther.mis_defined_SpecificHeatM6Var;
        mis_defined_SpecificHeatM7Var = rOther.mis_defined_SpecificHeatM7Var;
        mis_defined_ReactionVar = rOther.mis_defined_ReactionVar;
        mIsDefinedReactionGradientVar = rOther.mIsDefinedReactionGradientVar;

        return *this;
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ConvectionDiffusionSettings #" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConvectionDiffusionSettings #";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const Variable<double>* mpDensityVar = nullptr;
    //多相
    const Variable<double>* mpDensityM1Var = nullptr;
    const Variable<double>* mpDensityM2Var = nullptr;
    const Variable<double>* mpDensityM3Var = nullptr;
    const Variable<double>* mpDensityM4Var = nullptr;
    const Variable<double>* mpDensityM5Var = nullptr;
    const Variable<double>* mpDensityM6Var = nullptr;
    const Variable<double>* mpDensityM7Var = nullptr;
    const Variable<double>* mpMassFractionM1Var = nullptr;
    const Variable<double>* mpMassFractionM2Var = nullptr;
    const Variable<double>* mpMassFractionM3Var = nullptr;
    const Variable<double>* mpMassFractionM4Var = nullptr;
    const Variable<double>* mpMassFractionM5Var = nullptr;
    const Variable<double>* mpMassFractionM6Var = nullptr;
    const Variable<double>* mpMassFractionM7Var = nullptr;
    
    const Variable<double>* mpDiffusionVar = nullptr;
    const Variable<double>* mpDiffusionM1Var = nullptr;
    const Variable<double>* mpDiffusionM2Var = nullptr;
    const Variable<double>* mpDiffusionM3Var = nullptr;
    const Variable<double>* mpDiffusionM4Var = nullptr;
    const Variable<double>* mpDiffusionM5Var = nullptr;
    const Variable<double>* mpDiffusionM6Var = nullptr;
    const Variable<double>* mpDiffusionM7Var = nullptr;
    const Variable<double>* mpUnknownVar = nullptr;
    const Variable<double>* mpVolumeSourceVar = nullptr;
    const Variable<double>* mpAustenizeLatentVar= nullptr;
    const Variable<double>* mpEutectoidLatentVar= nullptr;
    const Variable<double>* mpBainiteLatentVar= nullptr;
    const Variable<double>* mpMartensiteLatentVar= nullptr;
    const Variable<double>* mpSurfaceSourceVar = nullptr;
    const Variable<double>* mpProjectionVar = nullptr;
    const Variable<array_1d<double,3> >* mpConvectionVar = nullptr;
    const Variable<array_1d<double,3> >* mpGradientVar = nullptr;
    const Variable<array_1d<double,3> >* mpMeshVelocityVar = nullptr;
    const Variable<double>* mpTransferCoefficientVar = nullptr;
    const Variable<array_1d<double,3> >* mpVelocityVar = nullptr;
    const Variable<double>* mpSpecificHeatVar = nullptr;
    const Variable<double>* mpSpecificHeatM1Var = nullptr;
    const Variable<double>* mpSpecificHeatM2Var = nullptr;
    const Variable<double>* mpSpecificHeatM3Var = nullptr;
    const Variable<double>* mpSpecificHeatM4Var = nullptr;
    const Variable<double>* mpSpecificHeatM5Var = nullptr;
    const Variable<double>* mpSpecificHeatM6Var = nullptr;
    const Variable<double>* mpSpecificHeatM7Var = nullptr;
    const Variable<double>* mpReactionVar = nullptr;
    const Variable<array_1d<double,3>>* mpReactionGradientVar = nullptr;
    bool mis_defined_DensityVar = false;
    //多相
    bool mis_defined_DensityM1Var = false;
    bool mis_defined_DensityM2Var = false;
    bool mis_defined_DensityM3Var = false;
    bool mis_defined_DensityM4Var = false;
    bool mis_defined_DensityM5Var = false;
    bool mis_defined_DensityM6Var = false;
    bool mis_defined_DensityM7Var = false;
    bool mis_defined_MassFractionM1Var = false;
    bool mis_defined_MassFractionM2Var = false;
    bool mis_defined_MassFractionM3Var = false;
    bool mis_defined_MassFractionM4Var = false;
    bool mis_defined_MassFractionM5Var = false;
    bool mis_defined_MassFractionM6Var = false;
    bool mis_defined_MassFractionM7Var = false;
    
    bool mis_defined_DiffusionVar = false;
    bool mis_defined_DiffusionM1Var = false;
    bool mis_defined_DiffusionM2Var = false;
    bool mis_defined_DiffusionM3Var = false;
    bool mis_defined_DiffusionM4Var = false;
    bool mis_defined_DiffusionM5Var = false;
    bool mis_defined_DiffusionM6Var = false;
    bool mis_defined_DiffusionM7Var = false;
    bool mis_defined_UnknownVar = false;
    bool mis_defined_VolumeSourceVar = false;
    bool mis_defined_AustenizeLatentVar = false;
    bool mis_defined_EutectoidLatentVar = false;
    bool mis_defined_BainiteLatentVar = false;
    bool mis_defined_MartensiteLatentVar = false;
    bool mis_defined_SurfaceSourceVar = false;
    bool mis_defined_ProjectionVar = false;
    bool mis_defined_ConvectionVar = false;
    bool mis_defined_GradientVar = false;
    bool mis_defined_MeshVelocityVar = false;
    bool mis_defined_TransferCoefficientVar = false;
    bool mis_defined_VelocityVar = false;
    bool mis_defined_SpecificHeatVar = false;
    bool mis_defined_SpecificHeatM1Var = false;
    bool mis_defined_SpecificHeatM2Var = false;
    bool mis_defined_SpecificHeatM3Var = false;
    bool mis_defined_SpecificHeatM4Var = false;
    bool mis_defined_SpecificHeatM5Var = false;
    bool mis_defined_SpecificHeatM6Var = false;
    bool mis_defined_SpecificHeatM7Var = false;
    bool mis_defined_ReactionVar = false;
    bool mIsDefinedReactionGradientVar = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;


    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{





    ///@}

}; // Class ConvectionDiffusionSettings

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ConvectionDiffusionSettings& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ConvectionDiffusionSettings& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(ConvectionDiffusionSettings::Pointer, CONVECTION_DIFFUSION_SETTINGS)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

}  // namespace Kratos.

#endif // KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED  defined
