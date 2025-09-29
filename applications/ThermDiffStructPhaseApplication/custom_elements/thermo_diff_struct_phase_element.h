// KRATOS MULTIPHYSICS COUPLING ELEMENT
// 
// Author: whf
// 
// Thermo-Chemo-Mechanical-PhaseChange Coupled Element

#pragma once

// System includes
#include <vector>
#include <array>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "../phase_transformation/phase_transformation_law.h"
#include "../phase_transformation/creep/creep_model.h"

namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(THERM_DIFF_STRUCT_PHASE_APPLICATION) ThermoDiffStructPhaseElement
    : public BaseSolidElement
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Base type definition
    typedef BaseSolidElement BaseType;
    
    /// Index type definition
    typedef std::size_t IndexType;
    
    /// Size type definition
    typedef std::size_t SizeType;

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;

    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// This is the definition of the node.
    typedef Node NodeType;
    
    /// Counted pointer definition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermoDiffStructPhaseElement);

    // Static constants
    static constexpr SizeType NumPhases = 7;
    static constexpr SizeType NumTransformations = 4;
    // static constexpr SizeType DofsPerNode = TDim + 2; // displacement + temperature + concentration
    // 定义每个节点的自由度数量
    static constexpr SizeType DofsPerNode3D = 5; // 3个位移 + 温度 + 浓度
    static constexpr SizeType DofsPerNode2D = 4; // 2个位移 + 温度 + 浓
    static constexpr SizeType TotalDofs2D = TNumNodes * 4;
    static constexpr SizeType TotalDofs3D = TNumNodes * 5;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ThermoDiffStructPhaseElement();
    /// Constructor with Id and geometry
    ThermoDiffStructPhaseElement(IndexType NewId, GeometryType::Pointer pGeometry);
    /// Constructor with Id, geometry and properties
    ThermoDiffStructPhaseElement(IndexType NewId, 
                                           GeometryType::Pointer pGeometry, 
                                           PropertiesType::Pointer pProperties);
    
    // copy constructor
    ThermoDiffStructPhaseElement(ThermoDiffStructPhaseElement const& rOther)
        : BaseType(rOther)
    {};

    /// Destructor
    virtual ~ThermoDiffStructPhaseElement() override = default;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, 
                           NodesArrayType const& ThisNodes, 
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ThermoDiffStructPhaseElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ThermoDiffStructPhaseElement>(
            NewId, pGeom, pProperties);
    }

    Element::Pointer Clone(IndexType NewId, 
                          NodesArrayType const& rThisNodes) const override;

    void EquationIdVector(EquationIdVectorType& rResult, 
                         const ProcessInfo& rCurrentProcessInfo) const override;
    
    void GetDofList(DofsVectorType& ElementalDofList, 
                   const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                             VectorType& rRightHandSideVector,
                             const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{
    
    // 返回单元信息
    std::string Info() const override
    {
        return "ThermoDiffStructPhaseElement";
    }

    // 返回单元名称字符串+节点和维度后缀
    std::string GetElementTypeAsString() const override
    {
        std::string BaseElementName = Element::GetElementTypeAsString();
        return "ThermoDiffStructPhase" + BaseElementName;
    }

    // 打印单元信息
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << TDim << "D" << TNumNodes << "N #" << " ID: " this->Id() <<
        "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

protected:
    ///@name Protected Structures
    ///@{
    
    struct ConvectionVariables
    {
        // Time integration parameters
        double theta;
        double dt;
        double dt_inv;
        double lumping_factor;
        
        // Current and previous step values
        array_1d<double, TNumNodes> temperature_current;
        array_1d<double, TNumNodes> temperature_previous;
        array_1d<double, TNumNodes> concentration_current;
        array_1d<double, TNumNodes> concentration_previous;
        
        // Phase fractions (7 phases for each node)
        // 0: Austenite, 1: Ferrite, 2: Cementite, 3: Pearlite, 4: Upper Bainite, 5: Lower Bainite, 6: Martensite
        array_1d<std::array<double, NumPhases>, TNumNodes> mass_fractions_current;
        array_1d<std::array<double, NumPhases>, TNumNodes> mass_fractions_previous;
        
        // Material properties (temperature and phase mass fraction dependent)
        double effective_density;
        double effective_specific_heat;
        double effective_conductivity;
        double effective_diffusivity;
        double effective_young_modulus;
        double effective_poisson_ratio;
        
        // Phase transformation latent
        array_1d<double, TNumNodes> phase_trans_latent_heat;
        array_1d<double, TNumNodes> austenization_latent;
        array_1d<double, TNumNodes> eutectoid_latent;
        array_1d<double, TNumNodes> bainitic_latent;
        array_1d<double, TNumNodes> martensitic_latent;


        //array_1d<double, TNumNodes> phase_transformation_strain[6]; // Voigt notation
        
        // Source terms
        array_1d<double, TNumNodes> thermal_source;
        array_1d<double, TNumNodes> diff_source;

        //default constructor
        ConvectionVariables(const SizeType StrainSize, const SizeType Dimension, const SizeType NumNodes)
        {
            theta = 0;
            dt = 0.0;
            dt_inv = 0.0;
            lumping_factor = 1.00 / double(NumNodes);
            effective_density = 0.0;
            effective_specific_heat = 0.0;
            effective_conductivity = 0.0;
            effective_diffusivity = 0.0;
            effective_young_modulus = 0.0;
            effective_poisson_ratio = 0.0;

            temperature_current = ZeroVector(NumNodes);
            temperature_previous = ZeroVector(NumNodes);
            concentration_current = ZeroVector(NumNodes);
            concentration_previous = ZeroVector(NumNodes);

            for (IndexType i = 0; i < NumNodes; ++i)
            {
                mass_fractions_current[i].fill(0.0);
                mass_fractions_previous[i].fill(0.0);
            }
            thermal_source = ZeroVector(NumNodes);
            diff_source = ZeroVector(NumNodes);

            phase_trans_latent_heat = ZeroVector(NumNodes);
            austenization_latent = ZeroVector(NumNodes);
            eutectoid_latent = ZeroVector(NumNodes);
            bainitic_latent = ZeroVector(NumNodes);
            martensitic_latent = ZeroVector(NumNodes);

        }
    };

    struct MechinicalVariables
    {
        Vector  N;
        Matrix  B;
        double  detF;
        Matrix  F;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;
        Vector Displacements_current;
        Vector Displacements_previous;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param Dimension The problem dimension: 2D or 3D
         * @param NumberOfNodes The number of nodes in the element
         */
        MechinicalVariables(
            const SizeType StrainSize,
            const SizeType Dimension,
            const SizeType NumberOfNodes
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
            Displacements_current = ZeroVector(Dimension * NumberOfNodes);
            Displacements_previous = ZeroVector(Dimension * NumberOfNodes);
        }
    };

    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const SizeType StrainSize)
        {
            if (StrainVector.size() != StrainSize)
                StrainVector.resize(StrainSize);

            if (StressVector.size() != StrainSize)
                StressVector.resize(StrainSize);

            if (D.size1() != StrainSize || D.size2() != StrainSize)
                D.resize(StrainSize, StrainSize);

            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    ///@}
    ///@name Protected Operations
    ///@{

    ThermoDiffStructPhaseElement() : BaseSolidElement()
    {
    }

    /**
     * @brief This method returns if the element provides the strain
     */
    bool UseElementProvidedStrain() const override;

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) override;

    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        ) override;

    /**
     * Calculation of the Deformation Matrix B
     * @param rB The deformation matrix
     * @param rDN_DX The derivatives of the shape functions
     * @param IntegrationPoints The array containing the integration points
     * @param PointNumber The integration point considered
     */
    virtual void CalculateB(
        Matrix& rB,
        const Matrix& rDN_DX,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const IndexType PointNumber
        ) const;

    /**
     * @brief Calculation of the equivalent deformation gradient
     * @param rF The deformation gradient F
     * @param StrainVector The strain tensor (Voigt notation)
     */
    virtual void ComputeEquivalentF(
        Matrix& rF,
        const Vector& StrainVector
        ) const;

    /**
     * @brief Initialize element variables
     */
    void InitializeCoupledElement(CoupledElementVariables& rVariables,
                                 const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Get nodal values for all fields
     */
    void GetAllNodalValues(CoupledElementVariables& rVariables,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate geometry derivatives and integration weights
     */
    void CalculateGeometryData(CoupledElementVariables& rVariables,
                              const IndexType IntegrationPoint);

    /**
     * @brief Update material properties based on phase composition
     */
    void UpdateMaterialProperties(CoupledElementVariables& rVariables,
                                 const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate phase transformations based on T and C
     */
    void CalculatePhaseTransformation(
        const CoupledElementVariables& rVariables,
        const IndexType NodeIndex,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate thermal contribution to system matrices
     */
    void CalculateThermalContribution(MatrixType& rLeftHandSideMatrix,
                                     VectorType& rRightHandSideVector,
                                     const CoupledElementVariables& rVariables,
                                     const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate diffusion contribution to system matrices
     */
    void CalculateDiffusionContribution(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       const CoupledElementVariables& rVariables,
                                       const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate mechanical contribution to system matrices
     */
    void CalculateMechanicalContribution(MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        const CoupledElementVariables& rVariables,
                                        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Add coupling terms between different fields
     */
    void AddCouplingTerms(MatrixType& rLeftHandSideMatrix,
                         VectorType& rRightHandSideVector,
                         const CoupledElementVariables& rVariables,
                         const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Assembly DOF ordering: [u1,v1,w1,T1,C1, u2,v2,w2,T2,C2, ...]
     */
    IndexType GetDofIndex(IndexType NodeIndex, IndexType DofType) const
    {
        // DofType: 0,1,2 = displacement components, 3 = temperature, 4 = concentration
        return NodeIndex * DofsPerNode + DofType;
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Calculate austenite transformation rate
     */
    double CalculateAusteniteTransformationRate(double temperature, 
                                               double concentration,
                                               const std::array<double, NumPhases>& current_phases);

    /**
     * @brief Calculate eutectoid transformation rate
     */
    double CalculateEutectoidTransformationRate(double temperature,
                                              double concentration,
                                              const std::array<double, NumPhases>& current_phases);

    /**
     * @brief Calculate bainite transformation rate
     */
    double CalculateBainiteTransformationRate(double temperature,
                                            double concentration,
                                            const std::array<double, NumPhases>& current_phases);

    /**
     * @brief Calculate martensite transformation rate
     */
    double CalculateMartensiteTransformationRate(double temperature,
                                               double concentration,
                                               const std::array<double, NumPhases>& current_phases);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseSolidElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseSolidElement);
    }

    ///@}

}; // Class ThermoDiffStructPhaseElement

} // namespace Kratos