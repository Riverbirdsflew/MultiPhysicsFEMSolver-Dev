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

// System includes

// External includes

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "includes/kratos_components.h"

namespace Kratos {

void ConvectionDiffusionSettings::save(Serializer& rSerializer) const
{
    KRATOS_TRY

    // Save the is defined bool flag for each variable
    rSerializer.save("mis_defined_DensityVar",mis_defined_DensityVar);
    //多相
    rSerializer.save("mis_defined_DensityM1Var",mis_defined_DensityM1Var);
    rSerializer.save("mis_defined_DensityM2Var",mis_defined_DensityM2Var);
    rSerializer.save("mis_defined_DensityM3Var",mis_defined_DensityM3Var);
    rSerializer.save("mis_defined_DensityM4Var",mis_defined_DensityM4Var);
    rSerializer.save("mis_defined_DensityM5Var",mis_defined_DensityM5Var);
    rSerializer.save("mis_defined_DensityM6Var",mis_defined_DensityM6Var);
    rSerializer.save("mis_defined_DensityM7Var",mis_defined_DensityM7Var);
    rSerializer.save("mis_defined_MassFractionM1Var",mis_defined_MassFractionM1Var);
    rSerializer.save("mis_defined_MassFractionM2Var",mis_defined_MassFractionM2Var);
    rSerializer.save("mis_defined_MassFractionM3Var",mis_defined_MassFractionM3Var);
    rSerializer.save("mis_defined_MassFractionM4Var",mis_defined_MassFractionM4Var);
    rSerializer.save("mis_defined_MassFractionM5Var",mis_defined_MassFractionM5Var);
    rSerializer.save("mis_defined_MassFractionM6Var",mis_defined_MassFractionM6Var);
    rSerializer.save("mis_defined_MassFractionM7Var",mis_defined_MassFractionM7Var);
    
    rSerializer.save("mis_defined_DiffusionVar",mis_defined_DiffusionVar);
    rSerializer.save("mis_defined_DiffusionM1Var",mis_defined_DiffusionM1Var);
    rSerializer.save("mis_defined_DiffusionM2Var",mis_defined_DiffusionM2Var);
    rSerializer.save("mis_defined_DiffusionM3Var",mis_defined_DiffusionM3Var);
    rSerializer.save("mis_defined_DiffusionM4Var",mis_defined_DiffusionM4Var);
    rSerializer.save("mis_defined_DiffusionM5Var",mis_defined_DiffusionM5Var);
    rSerializer.save("mis_defined_DiffusionM6Var",mis_defined_DiffusionM6Var);
    rSerializer.save("mis_defined_DiffusionM7Var",mis_defined_DiffusionM7Var);
    rSerializer.save("mis_defined_UnknownVar",mis_defined_UnknownVar);
    rSerializer.save("mis_defined_VolumeSourceVar",mis_defined_VolumeSourceVar);
    rSerializer.save("mis_defined_AustenizeLatentVar",mis_defined_AustenizeLatentVar);
    rSerializer.save("mis_defined_EutectoidLatentVar",mis_defined_EutectoidLatentVar);
    rSerializer.save("mis_defined_BainiteLatentVar",mis_defined_BainiteLatentVar);
    rSerializer.save("mis_defined_MartensiteLatentVar",mis_defined_MartensiteLatentVar);
    rSerializer.save("mis_defined_SurfaceSourceVar",mis_defined_SurfaceSourceVar);
    rSerializer.save("mis_defined_ProjectionVar",mis_defined_ProjectionVar);
    rSerializer.save("mis_defined_ConvectionVar",mis_defined_ConvectionVar);
    rSerializer.save("mis_defined_GradientVar",mis_defined_GradientVar);
    rSerializer.save("mis_defined_MeshVelocityVar",mis_defined_MeshVelocityVar);
    rSerializer.save("mis_defined_TransferCoefficientVar",mis_defined_TransferCoefficientVar);
    rSerializer.save("mis_defined_VelocityVar",mis_defined_VelocityVar);
    rSerializer.save("mis_defined_SpecificHeatVar",mis_defined_SpecificHeatVar);
    rSerializer.save("mis_defined_SpecificHeatM1Var",mis_defined_SpecificHeatM1Var);
    rSerializer.save("mis_defined_SpecificHeatM2Var",mis_defined_SpecificHeatM2Var);
    rSerializer.save("mis_defined_SpecificHeatM3Var",mis_defined_SpecificHeatM3Var);
    rSerializer.save("mis_defined_SpecificHeatM4Var",mis_defined_SpecificHeatM4Var);
    rSerializer.save("mis_defined_SpecificHeatM5Var",mis_defined_SpecificHeatM5Var);
    rSerializer.save("mis_defined_SpecificHeatM6Var",mis_defined_SpecificHeatM6Var);
    rSerializer.save("mis_defined_SpecificHeatM7Var",mis_defined_SpecificHeatM7Var);
    rSerializer.save("mis_defined_ReactionVar",mis_defined_ReactionVar);
    rSerializer.save("mIsDefinedReactionGradientVar", mIsDefinedReactionGradientVar);

    // Save the variable names
    // Note that the variable class save method only saves the name of the variables
    if (mpDensityVar != nullptr && mis_defined_DensityVar) {
        rSerializer.save("DensityVarName",mpDensityVar);
    }
    //多相
    if (mpDensityM1Var != nullptr && mis_defined_DensityM1Var) {
        rSerializer.save("DensityM1VarName",mpDensityM1Var);
    }
    if (mpDensityM2Var != nullptr && mis_defined_DensityM2Var) {
        rSerializer.save("DensityM2VarName",mpDensityM2Var);
    }
    if (mpDensityM3Var != nullptr && mis_defined_DensityM3Var) {
        rSerializer.save("DensityM3VarName",mpDensityM3Var);
    }
    if (mpDensityM4Var != nullptr && mis_defined_DensityM4Var) {
        rSerializer.save("DensityM4VarName",mpDensityM4Var);
    }
    if (mpDensityM5Var != nullptr && mis_defined_DensityM5Var) {
        rSerializer.save("DensityM5VarName",mpDensityM5Var);
    }
    if (mpDensityM6Var != nullptr && mis_defined_DensityM6Var) {
        rSerializer.save("DensityM6VarName",mpDensityM6Var);
    }
    if (mpDensityM7Var != nullptr && mis_defined_DensityM7Var) {
        rSerializer.save("DensityM7VarName",mpDensityM7Var);
    }
    if (mpMassFractionM1Var != nullptr && mis_defined_MassFractionM1Var) {
        rSerializer.save("MassFractionM1VarName",mpMassFractionM1Var);
    }
    if (mpMassFractionM2Var != nullptr && mis_defined_MassFractionM2Var) {
        rSerializer.save("MassFractionM2VarName",mpMassFractionM2Var);
    }
    if (mpMassFractionM3Var != nullptr && mis_defined_MassFractionM3Var) {
        rSerializer.save("MassFractionM3VarName",mpMassFractionM3Var);
    }
    if (mpMassFractionM4Var != nullptr && mis_defined_MassFractionM4Var) {
        rSerializer.save("MassFractionM4VarName",mpMassFractionM4Var);
    }
    if (mpMassFractionM5Var != nullptr && mis_defined_MassFractionM5Var) {
        rSerializer.save("MassFractionM5VarName",mpMassFractionM5Var);
    }
    if (mpMassFractionM6Var != nullptr && mis_defined_MassFractionM6Var) {
        rSerializer.save("MassFractionM6VarName",mpMassFractionM6Var);
    }
    if (mpMassFractionM7Var != nullptr && mis_defined_MassFractionM7Var) {
        rSerializer.save("MassFractionM7VarName",mpMassFractionM7Var);
    }
    
    if (mpDiffusionVar != nullptr && mis_defined_DiffusionVar) {
        rSerializer.save("DiffusionVarName",mpDiffusionVar);
    }
    if (mpDiffusionM1Var != nullptr && mis_defined_DiffusionM1Var) {
        rSerializer.save("DiffusionM1VarName",mpDiffusionM1Var);
    }
    if (mpDiffusionM2Var != nullptr && mis_defined_DiffusionM2Var) {
        rSerializer.save("DiffusionM2VarName",mpDiffusionM2Var);
    }
    if (mpDiffusionM3Var != nullptr && mis_defined_DiffusionM3Var) {
        rSerializer.save("DiffusionM3VarName",mpDiffusionM3Var);
    }
    if (mpDiffusionM4Var != nullptr && mis_defined_DiffusionM4Var) {
        rSerializer.save("DiffusionM4VarName",mpDiffusionM4Var);
    }
    if (mpDiffusionM5Var != nullptr && mis_defined_DiffusionM5Var) {
        rSerializer.save("DiffusionM5VarName",mpDiffusionM5Var);
    }
    if (mpDiffusionM6Var != nullptr && mis_defined_DiffusionM6Var) {
        rSerializer.save("DiffusionM6VarName",mpDiffusionM6Var);
    }
    if (mpDiffusionM7Var != nullptr && mis_defined_DiffusionM7Var) {
        rSerializer.save("DiffusionM7VarName",mpDiffusionM7Var);
    }
    if (mpUnknownVar != nullptr && mis_defined_UnknownVar) {
        rSerializer.save("UnknownVarName",mpUnknownVar);
    }
    if (mpVolumeSourceVar != nullptr && mis_defined_VolumeSourceVar) {
        rSerializer.save("VolumeSourceVarName",mpVolumeSourceVar);
    }
    if (mpAustenizeLatentVar != nullptr && mis_defined_AustenizeLatentVar) {
        rSerializer.save("AustenizeLatentVarName",mpAustenizeLatentVar);
    }
    if (mpEutectoidLatentVar != nullptr && mis_defined_EutectoidLatentVar) {
        rSerializer.save("EutectoidLatentVarName",mpEutectoidLatentVar);
    }
    if (mpBainiteLatentVar != nullptr && mis_defined_BainiteLatentVar) {
        rSerializer.save("BainiteLatentVarName",mpBainiteLatentVar);
    }
    if (mpMartensiteLatentVar != nullptr && mis_defined_MartensiteLatentVar) {
        rSerializer.save("MartensiteLatentVarName",mpMartensiteLatentVar);
    }
    if (mpSurfaceSourceVar != nullptr && mis_defined_SurfaceSourceVar) {
        rSerializer.save("SurfaceSourceVarName",mpSurfaceSourceVar);
    }
    if (mpProjectionVar != nullptr && mis_defined_ProjectionVar) {
        rSerializer.save("ProjectionVarName",mpProjectionVar);
    }
    if (mpConvectionVar != nullptr && mis_defined_ConvectionVar) {
        rSerializer.save("ConvectionVarName",mpConvectionVar);
    }
    if (mpGradientVar != nullptr && mis_defined_GradientVar) {
        rSerializer.save("GradientVarName",mpGradientVar);
    }
    if (mpMeshVelocityVar != nullptr && mis_defined_MeshVelocityVar) {
        rSerializer.save("MeshVelocityVarName",mpMeshVelocityVar);
    }
    if (mpTransferCoefficientVar != nullptr && mis_defined_TransferCoefficientVar) {
        rSerializer.save("TransferCoefficientVarName",mpTransferCoefficientVar);
    }
    if (mpVelocityVar != nullptr && mis_defined_VelocityVar) {
        rSerializer.save("VelocityVarName",mpVelocityVar);
    }
    if (mpSpecificHeatVar != nullptr && mis_defined_SpecificHeatVar) {
        rSerializer.save("SpecificHeatVarName",mpSpecificHeatVar);
    }
    if (mpSpecificHeatM1Var != nullptr && mis_defined_SpecificHeatM1Var) {
        rSerializer.save("SpecificHeatM1VarName",mpSpecificHeatM1Var);
    }
    if (mpSpecificHeatM2Var != nullptr && mis_defined_SpecificHeatM2Var) {
        rSerializer.save("SpecificHeatM2VarName",mpSpecificHeatM2Var);
    }
    if (mpSpecificHeatM3Var != nullptr && mis_defined_SpecificHeatM3Var) {
        rSerializer.save("SpecificHeatM3VarName",mpSpecificHeatM3Var);
    }
    if (mpSpecificHeatM4Var != nullptr && mis_defined_SpecificHeatM4Var) {
        rSerializer.save("SpecificHeatM4VarName",mpSpecificHeatM4Var);
    }
    if (mpSpecificHeatM5Var != nullptr && mis_defined_SpecificHeatM5Var) {
        rSerializer.save("SpecificHeatM5VarName",mpSpecificHeatM5Var);
    }
    if (mpSpecificHeatM6Var != nullptr && mis_defined_SpecificHeatM6Var) {
        rSerializer.save("SpecificHeatM6VarName",mpSpecificHeatM6Var);
    }
    if (mpSpecificHeatM7Var != nullptr && mis_defined_SpecificHeatM7Var) {
        rSerializer.save("SpecificHeatM7VarName",mpSpecificHeatM7Var);
    }
    if (mpReactionVar != nullptr && mis_defined_ReactionVar) {
        rSerializer.save("ReactionVarName",mpReactionVar);
    }
    if (mpReactionGradientVar != nullptr && mIsDefinedReactionGradientVar) {
        rSerializer.save("ReactionGradientVarName",mpReactionGradientVar);
    }

    KRATOS_CATCH("")
}

void ConvectionDiffusionSettings::load(Serializer& rSerializer)
{
    KRATOS_TRY

    // Load the is defined bool flags for each variable
    rSerializer.load("mis_defined_DensityVar",mis_defined_DensityVar);
    //多相
    rSerializer.load("mis_defined_DensityM1Var",mis_defined_DensityM1Var);
    rSerializer.load("mis_defined_DensityM2Var",mis_defined_DensityM2Var);
    rSerializer.load("mis_defined_DensityM3Var",mis_defined_DensityM3Var);
    rSerializer.load("mis_defined_DensityM4Var",mis_defined_DensityM4Var);
    rSerializer.load("mis_defined_DensityM5Var",mis_defined_DensityM5Var);
    rSerializer.load("mis_defined_DensityM6Var",mis_defined_DensityM6Var);
    rSerializer.load("mis_defined_DensityM7Var",mis_defined_DensityM7Var);
    rSerializer.load("mis_defined_MassFractionM1Var",mis_defined_MassFractionM1Var);
    rSerializer.load("mis_defined_MassFractionM2Var",mis_defined_MassFractionM2Var);
    rSerializer.load("mis_defined_MassFractionM3Var",mis_defined_MassFractionM3Var);
    rSerializer.load("mis_defined_MassFractionM4Var",mis_defined_MassFractionM4Var);
    rSerializer.load("mis_defined_MassFractionM5Var",mis_defined_MassFractionM5Var);
    rSerializer.load("mis_defined_MassFractionM6Var",mis_defined_MassFractionM6Var);
    rSerializer.load("mis_defined_MassFractionM7Var",mis_defined_MassFractionM7Var);
    
    rSerializer.load("mis_defined_DiffusionVar",mis_defined_DiffusionVar);
    rSerializer.load("mis_defined_DiffusionM1Var",mis_defined_DiffusionM1Var);
    rSerializer.load("mis_defined_DiffusionM2Var",mis_defined_DiffusionM2Var);
    rSerializer.load("mis_defined_DiffusionM3Var",mis_defined_DiffusionM3Var);
    rSerializer.load("mis_defined_DiffusionM4Var",mis_defined_DiffusionM4Var);
    rSerializer.load("mis_defined_DiffusionM5Var",mis_defined_DiffusionM5Var);
    rSerializer.load("mis_defined_DiffusionM6Var",mis_defined_DiffusionM6Var);
    rSerializer.load("mis_defined_DiffusionM7Var",mis_defined_DiffusionM7Var);
    rSerializer.load("mis_defined_UnknownVar",mis_defined_UnknownVar);
    rSerializer.load("mis_defined_VolumeSourceVar",mis_defined_VolumeSourceVar);
    rSerializer.load("mis_defined_AustenizeLatentVar",mis_defined_AustenizeLatentVar);
    rSerializer.load("mis_defined_EutectoidLatentVar",mis_defined_EutectoidLatentVar);
    rSerializer.load("mis_defined_BainiteLatentVar",mis_defined_BainiteLatentVar);
    rSerializer.load("mis_defined_MartensiteLatentVar",mis_defined_MartensiteLatentVar);
    rSerializer.load("mis_defined_SurfaceSourceVar",mis_defined_SurfaceSourceVar);
    rSerializer.load("mis_defined_ProjectionVar",mis_defined_ProjectionVar);
    rSerializer.load("mis_defined_ConvectionVar",mis_defined_ConvectionVar);
    rSerializer.load("mis_defined_GradientVar",mis_defined_GradientVar);
    rSerializer.load("mis_defined_MeshVelocityVar",mis_defined_MeshVelocityVar);
    rSerializer.load("mis_defined_TransferCoefficientVar",mis_defined_TransferCoefficientVar);
    rSerializer.load("mis_defined_VelocityVar",mis_defined_VelocityVar);
    rSerializer.load("mis_defined_SpecificHeatVar",mis_defined_SpecificHeatVar);
    rSerializer.load("mis_defined_SpecificHeatM1Var",mis_defined_SpecificHeatM1Var);
    rSerializer.load("mis_defined_SpecificHeatM2Var",mis_defined_SpecificHeatM2Var);
    rSerializer.load("mis_defined_SpecificHeatM3Var",mis_defined_SpecificHeatM3Var);
    rSerializer.load("mis_defined_SpecificHeatM4Var",mis_defined_SpecificHeatM4Var);
    rSerializer.load("mis_defined_SpecificHeatM5Var",mis_defined_SpecificHeatM5Var);
    rSerializer.load("mis_defined_SpecificHeatM6Var",mis_defined_SpecificHeatM6Var);
    rSerializer.load("mis_defined_SpecificHeatM7Var",mis_defined_SpecificHeatM7Var);
    rSerializer.load("mis_defined_ReactionVar",mis_defined_ReactionVar);
    rSerializer.load("mIsDefinedReactionGradientVar", mIsDefinedReactionGradientVar);

    // If the variables are defined, load their name
    // Note that only the name has been saved to retrieve the already existent variable from the KratosComponents
    if(mis_defined_DensityVar) {
        std::string density_var_name;
        rSerializer.load("DensityVarName", density_var_name);
        mpDensityVar = &(KratosComponents<Variable<double>>::Get(density_var_name));
    }
    //多相
    if(mis_defined_DensityM1Var) {
        std::string density_m1_var_name;
        rSerializer.load("DensityM1VarName", density_m1_var_name);
        mpDensityM1Var = &(KratosComponents<Variable<double>>::Get(density_m1_var_name));
    }
    if(mis_defined_DensityM2Var) {
        std::string density_m2_var_name;
        rSerializer.load("DensityM2VarName", density_m2_var_name);
        mpDensityM2Var = &(KratosComponents<Variable<double>>::Get(density_m2_var_name));
    }
    if(mis_defined_DensityM3Var) {
        std::string density_m3_var_name;
        rSerializer.load("DensityM3VarName", density_m3_var_name);
        mpDensityM3Var = &(KratosComponents<Variable<double>>::Get(density_m3_var_name));
    }
    if(mis_defined_DensityM4Var) {
        std::string density_m4_var_name;
        rSerializer.load("DensityM4VarName", density_m4_var_name);
        mpDensityM4Var = &(KratosComponents<Variable<double>>::Get(density_m4_var_name));
    }
    if(mis_defined_DensityM5Var) {
        std::string density_m5_var_name;
        rSerializer.load("DensityM5VarName", density_m5_var_name);
        mpDensityM5Var = &(KratosComponents<Variable<double>>::Get(density_m5_var_name));
    }
    if(mis_defined_DensityM6Var) {
        std::string density_m6_var_name;
        rSerializer.load("DensityM6VarName", density_m6_var_name);
        mpDensityM6Var = &(KratosComponents<Variable<double>>::Get(density_m6_var_name));
    }
    if(mis_defined_DensityM7Var) {
        std::string density_m7_var_name;
        rSerializer.load("DensityM7VarName", density_m7_var_name);
        mpDensityM7Var = &(KratosComponents<Variable<double>>::Get(density_m7_var_name));
    }
    if(mis_defined_MassFractionM1Var) {
        std::string mass_fraction_m1_var_name;
        rSerializer.load("MassFractionM1VarName", mass_fraction_m1_var_name);
        mpMassFractionM1Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m1_var_name));
    }
    if(mis_defined_MassFractionM2Var) {
        std::string mass_fraction_m2_var_name;
        rSerializer.load("MassFractionM2VarName", mass_fraction_m2_var_name);
        mpMassFractionM2Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m2_var_name));
    }
    if(mis_defined_MassFractionM3Var) {
        std::string mass_fraction_m3_var_name;
        rSerializer.load("MassFractionM3VarName", mass_fraction_m3_var_name);
        mpMassFractionM3Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m3_var_name));
    }
    if(mis_defined_MassFractionM4Var) {
        std::string mass_fraction_m4_var_name;
        rSerializer.load("MassFractionM4VarName", mass_fraction_m4_var_name);
        mpMassFractionM4Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m4_var_name));
    }
    if(mis_defined_MassFractionM5Var) {
        std::string mass_fraction_m5_var_name;
        rSerializer.load("MassFractionM5VarName", mass_fraction_m5_var_name);
        mpMassFractionM5Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m5_var_name));
    }
    if(mis_defined_MassFractionM6Var) {
        std::string mass_fraction_m6_var_name;
        rSerializer.load("MassFractionM6VarName", mass_fraction_m6_var_name);
        mpMassFractionM6Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m6_var_name));
    }
    if(mis_defined_MassFractionM7Var) {
        std::string mass_fraction_m7_var_name;
        rSerializer.load("MassFractionM7VarName", mass_fraction_m7_var_name);
        mpMassFractionM7Var = &(KratosComponents<Variable<double>>::Get(mass_fraction_m7_var_name));
    }//多相
    
    if(mis_defined_DiffusionVar) {
        std::string diffusion_var_name;
        rSerializer.load("DiffusionVarName", diffusion_var_name);
        mpDiffusionVar = &(KratosComponents<Variable<double>>::Get(diffusion_var_name));
    }
    if(mis_defined_DiffusionM1Var) {
        std::string diffusion_m1_var_name;
        rSerializer.load("DiffusionM1VarName", diffusion_m1_var_name);
        mpDiffusionM1Var = &(KratosComponents<Variable<double>>::Get(diffusion_m1_var_name));
    }
    if(mis_defined_DiffusionM2Var) {
        std::string diffusion_m2_var_name;
        rSerializer.load("DiffusionM2VarName", diffusion_m2_var_name);
        mpDiffusionM2Var = &(KratosComponents<Variable<double>>::Get(diffusion_m2_var_name));
    }
    if(mis_defined_DiffusionM3Var) {
        std::string diffusion_m3_var_name;
        rSerializer.load("DiffusionM3VarName", diffusion_m3_var_name);
        mpDiffusionM3Var = &(KratosComponents<Variable<double>>::Get(diffusion_m3_var_name));
    }
    if(mis_defined_DiffusionM4Var) {
        std::string diffusion_m4_var_name;
        rSerializer.load("DiffusionM4VarName", diffusion_m4_var_name);
        mpDiffusionM4Var = &(KratosComponents<Variable<double>>::Get(diffusion_m4_var_name));
    }
    if(mis_defined_DiffusionM5Var) {
        std::string diffusion_m5_var_name;
        rSerializer.load("DiffusionM5VarName", diffusion_m5_var_name);
        mpDiffusionM5Var = &(KratosComponents<Variable<double>>::Get(diffusion_m5_var_name));
    }
    if(mis_defined_DiffusionM6Var) {
        std::string diffusion_m6_var_name;
        rSerializer.load("DiffusionM6VarName", diffusion_m6_var_name);
        mpDiffusionM6Var = &(KratosComponents<Variable<double>>::Get(diffusion_m6_var_name));
    }
    if(mis_defined_DiffusionM7Var) {
        std::string diffusion_m7_var_name;
        rSerializer.load("DiffusionM7VarName", diffusion_m7_var_name);
        mpDiffusionM7Var = &(KratosComponents<Variable<double>>::Get(diffusion_m7_var_name));
    }
    
    if(mis_defined_UnknownVar) {
        std::string unknown_var_name;
        rSerializer.load("UnknownVarName", unknown_var_name);
        mpUnknownVar = &(KratosComponents<Variable<double>>::Get(unknown_var_name));
    }
    if(mis_defined_VolumeSourceVar) {
        std::string volume_source_var_name;
        rSerializer.load("VolumeSourceVarName", volume_source_var_name);
        mpVolumeSourceVar = &(KratosComponents<Variable<double>>::Get(volume_source_var_name));
    }
    if(mis_defined_AustenizeLatentVar) {
        std::string austenize_latent_var_name;
        rSerializer.load("AustenizeLatentVarName", austenize_latent_var_name);
        mpAustenizeLatentVar = &(KratosComponents<Variable<double>>::Get(austenize_latent_var_name));
    }
    if(mis_defined_EutectoidLatentVar) {
        std::string eutectoid_latent_var_name;
        rSerializer.load("EutectoidLatentVarName", eutectoid_latent_var_name);
        mpEutectoidLatentVar = &(KratosComponents<Variable<double>>::Get(eutectoid_latent_var_name));
    }
    if(mis_defined_BainiteLatentVar) {
        std::string bainite_latent_var_name;
        rSerializer.load("BainiteLatentVarName", bainite_latent_var_name);
        mpBainiteLatentVar = &(KratosComponents<Variable<double>>::Get(bainite_latent_var_name));
    }
    if(mis_defined_MartensiteLatentVar) {
        std::string martensite_latent_var_name;
        rSerializer.load("MartensiteLatentVarName", martensite_latent_var_name);
        mpMartensiteLatentVar = &(KratosComponents<Variable<double>>::Get(martensite_latent_var_name));
    }
    if(mis_defined_SurfaceSourceVar) {
        std::string surface_source_var_name;
        rSerializer.load("SurfaceSourceVarName", surface_source_var_name);
        mpSurfaceSourceVar = &(KratosComponents<Variable<double>>::Get(surface_source_var_name));
    }
    if(mis_defined_ProjectionVar) {
        std::string projection_var_name;
        rSerializer.load("ProjectionVarName", projection_var_name);
        mpProjectionVar = &(KratosComponents<Variable<double>>::Get(projection_var_name));
    }
    if(mis_defined_ConvectionVar) {
        std::string convection_var_name;
        rSerializer.load("ConvectionVarName", convection_var_name);
        mpConvectionVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(convection_var_name));
    }
    if(mis_defined_GradientVar) {
        std::string gradient_var_name;
        rSerializer.load("GradientVarName", gradient_var_name);
        mpGradientVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(gradient_var_name));
    }
    if(mis_defined_MeshVelocityVar) {
        std::string mesh_velocity_var;
        rSerializer.load("MeshVelocityVarName", mesh_velocity_var);
        mpMeshVelocityVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(mesh_velocity_var));
    }
    if(mis_defined_TransferCoefficientVar) {
        std::string transfer_coefficient_var_name;
        rSerializer.load("TransferCoefficientVarName", transfer_coefficient_var_name);
        mpTransferCoefficientVar = &(KratosComponents<Variable<double>>::Get(transfer_coefficient_var_name));
    }
    if(mis_defined_VelocityVar) {
        std::string velocity_var_name;
        rSerializer.load("VelocityVarName", velocity_var_name);
        mpVelocityVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(velocity_var_name));
    }
    if(mis_defined_SpecificHeatVar) {
        std::string specific_heat_var_name;
        rSerializer.load("SpecificHeatVarName", specific_heat_var_name);
        mpSpecificHeatVar = &(KratosComponents<Variable<double>>::Get(specific_heat_var_name));
    }
    if(mis_defined_SpecificHeatM1Var) {
        std::string specific_heat_m1_var_name;
        rSerializer.load("SpecificHeatM1VarName", specific_heat_m1_var_name);
        mpSpecificHeatM1Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m1_var_name));
    }
    if(mis_defined_SpecificHeatM2Var) {
        std::string specific_heat_m2_var_name;
        rSerializer.load("SpecificHeatM2VarName", specific_heat_m2_var_name);
        mpSpecificHeatM2Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m2_var_name));
    }
    if(mis_defined_SpecificHeatM3Var) {
        std::string specific_heat_m3_var_name;
        rSerializer.load("SpecificHeatM3VarName", specific_heat_m3_var_name);
        mpSpecificHeatM3Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m3_var_name));
    }
    if(mis_defined_SpecificHeatM4Var) {
        std::string specific_heat_m4_var_name;
        rSerializer.load("SpecificHeatM4VarName", specific_heat_m4_var_name);
        mpSpecificHeatM4Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m4_var_name));
    }
    if(mis_defined_SpecificHeatM5Var) {
        std::string specific_heat_m5_var_name;
        rSerializer.load("SpecificHeatM5VarName", specific_heat_m5_var_name);
        mpSpecificHeatM5Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m5_var_name));
    }
    if(mis_defined_SpecificHeatM6Var) {
        std::string specific_heat_m6_var_name;
        rSerializer.load("SpecificHeatM6VarName", specific_heat_m6_var_name);
        mpSpecificHeatM6Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m6_var_name));
    }
    if(mis_defined_SpecificHeatM7Var) {
        std::string specific_heat_m7_var_name;
        rSerializer.load("SpecificHeatM7VarName", specific_heat_m7_var_name);
        mpSpecificHeatM7Var = &(KratosComponents<Variable<double>>::Get(specific_heat_m7_var_name));
    }
    
    if(mis_defined_ReactionVar) {
        std::string reaction_var_name;
        rSerializer.load("ReactionVarName", reaction_var_name);
        mpReactionVar = &(KratosComponents<Variable<double>>::Get(reaction_var_name));
    }
    if(mIsDefinedReactionGradientVar) {
        std::string reaction_gradient_var_name;
        rSerializer.load("ReactionGradientVarName", reaction_gradient_var_name);
        mpReactionGradientVar = &(KratosComponents<Variable<array_1d<double,3>>>::Get(reaction_gradient_var_name));
    }

    KRATOS_CATCH("")
}

}  // namespace Kratos.
