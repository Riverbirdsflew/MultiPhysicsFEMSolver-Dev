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

#if !defined(KRATOS_CONCENTRATION_DIFFUSION_SETTINGS_INCLUDED )
#define  KRATOS_CONCENTRATION_DIFFUSION_SETTINGS_INCLUDED

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

/// Concentration diffusion settings. This class contains information to be used by the Concentration diffusion solver, all the variables that will be needed by the solver.
/** All the variables needed by any Concentration diffusion problem are included here. However, many problems do not require such a large amount of variables.
 * For those cases, variables that are not defined will be set to either zero or one (depending on the variable)
 * For this purpose, there are flags to ask whether a variable has been defined or not. For each variable, there are three main functions.
 * SetVariableforXuse: we assign the variable for this use. When doing that we also set the flag that now this variable has been defined.
 * GetVariableforXuse: we return the variable for this use.
 * IsDefinedVariableforXuse: tells whether that variable has been defined or not.
*/
class KRATOS_API(KRATOS_CORE) ConcentrationDiffusionSettings
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConcentrationDiffusionSettings
    KRATOS_CLASS_POINTER_DEFINITION(ConcentrationDiffusionSettings);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConcentrationDiffusionSettings() = default;

    ConcentrationDiffusionSettings(const ConcentrationDiffusionSettings& rOther):
        //mpDensityVar(rOther.mpDensityVar),//密度
        mpMassFractionM1Var(rOther.mpMassFractionM1Var),//多相质量分数7
        mpMassFractionM2Var(rOther.mpMassFractionM2Var),
        mpMassFractionM3Var(rOther.mpMassFractionM3Var),
        mpMassFractionM4Var(rOther.mpMassFractionM4Var),
        mpMassFractionM5Var(rOther.mpMassFractionM5Var),
        mpMassFractionM6Var(rOther.mpMassFractionM6Var),
        mpMassFractionM7Var(rOther.mpMassFractionM7Var),
        mpDiffusionVar(rOther.mpDiffusionVar),
        mpDiffusionM1Var(rOther.mpDiffusionM1Var),//多相扩散系数7
        mpDiffusionM2Var(rOther.mpDiffusionM2Var),
        mpDiffusionM3Var(rOther.mpDiffusionM3Var),
        mpDiffusionM4Var(rOther.mpDiffusionM4Var),
        mpDiffusionM5Var(rOther.mpDiffusionM5Var),
        mpDiffusionM6Var(rOther.mpDiffusionM6Var),
        mpDiffusionM7Var(rOther.mpDiffusionM7Var),
        mpUnknownVar(rOther.mpUnknownVar),
        mpPhaseFractionVar(rOther.mpPhaseFractionVar),
        mpPhaseNameVar(rOther.mpPhaseNameVar),
        mpVolumeSourceVar(rOther.mpVolumeSourceVar),
        mpSurfaceSourceVar(rOther.mpSurfaceSourceVar),
        mpProjectionVar(rOther. mpProjectionVar),
        mpConcentrationVar(rOther.mpConcentrationVar),
        mpGradientVar(rOther.mpGradientVar),
        mpMeshVelocityVar(rOther.mpMeshVelocityVar),
        mpTransferCoefficientVar(rOther.mpTransferCoefficientVar),
        mpVelocityVar(rOther.mpVelocityVar),
        //mpSpecificHeatVar(rOther.mpSpecificHeatVar),//比热
        mpReactionVar(rOther.mpReactionVar),
        mpReactionGradientVar(rOther.mpReactionGradientVar),
        //mis_defined_DensityVar(rOther.mis_defined_DensityVar),//密度
		mis_defined_MassFractionM1Var(rOther.mis_defined_MassFractionM1Var),//多相扩散系数
        mis_defined_MassFractionM2Var(rOther.mis_defined_MassFractionM2Var),
        mis_defined_MassFractionM3Var(rOther.mis_defined_MassFractionM3Var),
        mis_defined_MassFractionM4Var(rOther.mis_defined_MassFractionM4Var),
        mis_defined_MassFractionM5Var(rOther.mis_defined_MassFractionM5Var),
        mis_defined_MassFractionM6Var(rOther.mis_defined_MassFractionM6Var),
        mis_defined_MassFractionM7Var(rOther.mis_defined_MassFractionM7Var),
        mis_defined_DiffusionVar(rOther.mis_defined_DiffusionVar),
        mis_defined_DiffusionM1Var(rOther.mis_defined_DiffusionM1Var),//多相扩散系数
        mis_defined_DiffusionM2Var(rOther.mis_defined_DiffusionM2Var),
        mis_defined_DiffusionM3Var(rOther.mis_defined_DiffusionM3Var),
        mis_defined_DiffusionM4Var(rOther.mis_defined_DiffusionM4Var),
        mis_defined_DiffusionM5Var(rOther.mis_defined_DiffusionM5Var),
        mis_defined_DiffusionM6Var(rOther.mis_defined_DiffusionM6Var),
        mis_defined_DiffusionM7Var(rOther.mis_defined_DiffusionM7Var),
		mis_defined_UnknownVar(rOther.mis_defined_UnknownVar),
		mis_defined_PhaseFractionVar(rOther.mis_defined_PhaseFractionVar),
		mis_defined_PhaseNameVar(rOther.mis_defined_PhaseNameVar),
		mis_defined_VolumeSourceVar(rOther.mis_defined_VolumeSourceVar),
		mis_defined_SurfaceSourceVar(rOther.mis_defined_SurfaceSourceVar),
		mis_defined_ProjectionVar(rOther.mis_defined_ProjectionVar),
		mis_defined_ConcentrationVar(rOther.mis_defined_ConcentrationVar),
        mis_defined_GradientVar(rOther.mis_defined_GradientVar),
		mis_defined_MeshVelocityVar(rOther.mis_defined_MeshVelocityVar),
		mis_defined_TransferCoefficientVar(rOther.mis_defined_TransferCoefficientVar),
		mis_defined_VelocityVar(rOther.mis_defined_VelocityVar),
		//mis_defined_SpecificHeatVar(rOther.mis_defined_SpecificHeatVar),//比热
        mis_defined_ReactionVar(rOther.mis_defined_ReactionVar),
        mIsDefinedReactionGradientVar(rOther.mIsDefinedReactionGradientVar)
    {
    }

    /// Destructor.
    virtual ~ConcentrationDiffusionSettings() {};

    ///@}
    ///@name Operators
    ///@{
    //~ void SetDensityVariable(const Variable<double>& rvar)
    //~ {
        //~ mpDensityVar = &rvar;
        //~ mis_defined_DensityVar=true;
    //~ }
    //~ const Variable<double>& GetDensityVariable() const
    //~ {
        //~ return *mpDensityVar;
    //~ }
    //~ bool IsDefinedDensityVariable() const
    //~ {
		//~ return mpDensityVar != nullptr;
	//~ } //设置密度

    //多相材料质量分数
    //m1
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
    //多相材料扩散系数
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
    
    const Variable<double>& GetDiffusionVariable() const
    {
        return *mpDiffusionVar;
    }
    bool IsDefinedDiffusionVariable() const
    {
		return mpDiffusionVar != nullptr;
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
	
    void SetPhaseFractionVariable(const Variable<double>& rvar)
    {
        mpPhaseFractionVar = &rvar;
		mis_defined_PhaseFractionVar=true;
    }
    const Variable<double>& GetPhaseFractionVariable() const
    {
        return *mpPhaseFractionVar;
    }
    bool IsDefinedPhaseFractionVariable() const
    {
		return mpPhaseFractionVar != nullptr;
	}

    void SetPhaseNameVariable(const Variable<double>& rvar)
    {
        mpPhaseNameVar = &rvar;
		mis_defined_PhaseNameVar=true;
    }
    const Variable<double>& GetPhaseNameVariable() const
    {
        return *mpPhaseNameVar;
    }
    bool IsDefinedPhaseNameVariable() const
    {
		return mpPhaseNameVar != nullptr;
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

    void SetConcentrationVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpConcentrationVar = &rvar;
		mis_defined_ConcentrationVar=true;
    }
    const Variable<array_1d<double,3> >& GetConcentrationVariable() const
    {
        return *mpConcentrationVar;
    }
    bool IsDefinedConcentrationVariable() const
    {
		return mpConcentrationVar != nullptr;
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

    //~ void SetSpecificHeatVariable(const Variable<double>& rvar)
    //~ {
        //~ mpSpecificHeatVar = &rvar;
		//~ mis_defined_SpecificHeatVar=true;
    //~ }
    //~ const Variable<double>& GetSpecificHeatVariable() const
    //~ {
        //~ return *mpSpecificHeatVar;
    //~ }
    //~ bool IsDefinedSpecificHeatVariable() const
    //~ {
		//~ return mpSpecificHeatVar != nullptr;
	//~ } //设置比热

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
    ConcentrationDiffusionSettings& operator=(ConcentrationDiffusionSettings const& rOther)
    {
        //mpDensityVar = rOther.mpDensityVar;//密度
        mpMassFractionM1Var = rOther.mpMassFractionM1Var;
        mpMassFractionM2Var = rOther.mpMassFractionM2Var;
        mpMassFractionM3Var = rOther.mpMassFractionM3Var;
        mpMassFractionM4Var = rOther.mpMassFractionM4Var;
        mpMassFractionM5Var = rOther.mpMassFractionM5Var;
        mpMassFractionM6Var = rOther.mpMassFractionM6Var;
        mpMassFractionM7Var = rOther.mpMassFractionM7Var;//多相
        mpDiffusionVar = rOther.mpDiffusionVar;
        mpDiffusionM1Var = rOther.mpDiffusionM1Var;//扩散系数
        mpDiffusionM2Var = rOther.mpDiffusionM2Var;
        mpDiffusionM3Var = rOther.mpDiffusionM3Var;
        mpDiffusionM4Var = rOther.mpDiffusionM4Var;
        mpDiffusionM5Var = rOther.mpDiffusionM5Var;
        mpDiffusionM6Var = rOther.mpDiffusionM6Var;
        mpDiffusionM7Var = rOther.mpDiffusionM7Var;
        mpUnknownVar = rOther.mpUnknownVar;
        mpPhaseFractionVar = rOther.mpPhaseFractionVar;
        mpPhaseNameVar = rOther.mpPhaseNameVar;
        mpVolumeSourceVar = rOther.mpVolumeSourceVar;
        mpSurfaceSourceVar = rOther.mpSurfaceSourceVar;
        mpProjectionVar = rOther.mpProjectionVar;
        mpConcentrationVar = rOther.mpConcentrationVar;
        mpGradientVar = rOther.mpGradientVar;
        mpMeshVelocityVar = rOther.mpMeshVelocityVar;
        mpTransferCoefficientVar = rOther.mpTransferCoefficientVar;
        mpVelocityVar = rOther.mpVelocityVar;
		//mpSpecificHeatVar = rOther.mpSpecificHeatVar;//比热
        mpReactionVar = rOther.mpReactionVar;
        mpReactionGradientVar = rOther.mpReactionGradientVar;
        //now the is_defined
        //mis_defined_DensityVar = rOther.mis_defined_DensityVar;//密度
		mis_defined_MassFractionM1Var = rOther.mis_defined_MassFractionM1Var;
        mis_defined_MassFractionM2Var = rOther.mis_defined_MassFractionM2Var;
        mis_defined_MassFractionM3Var = rOther.mis_defined_MassFractionM3Var;
        mis_defined_MassFractionM4Var = rOther.mis_defined_MassFractionM4Var;
        mis_defined_MassFractionM5Var = rOther.mis_defined_MassFractionM5Var;
        mis_defined_MassFractionM6Var = rOther.mis_defined_MassFractionM6Var;
        mis_defined_MassFractionM7Var = rOther.mis_defined_MassFractionM7Var;//多相
        mis_defined_DiffusionVar = rOther.mis_defined_DiffusionVar;
        mis_defined_DiffusionM1Var = rOther.mis_defined_DiffusionM1Var;//扩散系数
        mis_defined_DiffusionM2Var = rOther.mis_defined_DiffusionM2Var;
        mis_defined_DiffusionM3Var = rOther.mis_defined_DiffusionM3Var;
        mis_defined_DiffusionM4Var = rOther.mis_defined_DiffusionM4Var;
        mis_defined_DiffusionM5Var = rOther.mis_defined_DiffusionM5Var;
        mis_defined_DiffusionM6Var = rOther.mis_defined_DiffusionM6Var;
        mis_defined_DiffusionM7Var = rOther.mis_defined_DiffusionM7Var;
		mis_defined_UnknownVar = rOther.mis_defined_UnknownVar;
		mis_defined_PhaseFractionVar = rOther.mis_defined_PhaseFractionVar;
		mis_defined_PhaseNameVar = rOther.mis_defined_PhaseNameVar;
		mis_defined_VolumeSourceVar = rOther.mis_defined_VolumeSourceVar;
		mis_defined_SurfaceSourceVar = rOther.mis_defined_SurfaceSourceVar;
		mis_defined_ProjectionVar = rOther.mis_defined_ProjectionVar;
        mis_defined_ConcentrationVar = rOther.mis_defined_ConcentrationVar;
        mis_defined_GradientVar = rOther.mis_defined_GradientVar;
		mis_defined_MeshVelocityVar = rOther.mis_defined_MeshVelocityVar;
		mis_defined_TransferCoefficientVar = rOther.mis_defined_TransferCoefficientVar;
		mis_defined_VelocityVar = rOther.mis_defined_VelocityVar;
		//mis_defined_SpecificHeatVar = rOther.mis_defined_SpecificHeatVar;//比热
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
        buffer << "ConcentrationDiffusionSettings #" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConcentrationDiffusionSettings #";
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

    //const Variable<double>* mpDensityVar = nullptr;//密度
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
    const Variable<double>* mpPhaseFractionVar = nullptr;
    const Variable<double>* mpPhaseNameVar = nullptr;
    const Variable<double>* mpVolumeSourceVar = nullptr;
    const Variable<double>* mpSurfaceSourceVar = nullptr;
    const Variable<double>* mpProjectionVar = nullptr;
    const Variable<array_1d<double,3> >* mpConcentrationVar = nullptr;
    const Variable<array_1d<double,3> >* mpGradientVar = nullptr;
    const Variable<array_1d<double,3> >* mpMeshVelocityVar = nullptr;
    const Variable<double>* mpTransferCoefficientVar = nullptr;
    const Variable<array_1d<double,3> >* mpVelocityVar = nullptr;
    //const Variable<double>* mpSpecificHeatVar = nullptr;//比热
    const Variable<double>* mpReactionVar = nullptr;
    const Variable<array_1d<double,3>>* mpReactionGradientVar = nullptr;
    //bool mis_defined_DensityVar = false;//密度
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
    bool mis_defined_PhaseFractionVar = false;
    bool mis_defined_PhaseNameVar = false;
    bool mis_defined_VolumeSourceVar = false;
    bool mis_defined_SurfaceSourceVar = false;
    bool mis_defined_ProjectionVar = false;
    bool mis_defined_ConcentrationVar = false;
    bool mis_defined_GradientVar = false;
    bool mis_defined_MeshVelocityVar = false;
    bool mis_defined_TransferCoefficientVar = false;
    bool mis_defined_VelocityVar = false;
    //bool mis_defined_SpecificHeatVar = false;//比热
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

}; // Class ConcentrationDiffusionSettings

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ConcentrationDiffusionSettings& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ConcentrationDiffusionSettings& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(ConcentrationDiffusionSettings::Pointer, CONCENTRATION_DIFFUSION_SETTINGS)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

}  // namespace Kratos.

#endif // KRATOS_CONCENTRATION_DIFFUSION_SETTINGS_INCLUDED  defined
