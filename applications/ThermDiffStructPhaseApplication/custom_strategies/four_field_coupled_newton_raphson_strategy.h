#if !defined(KRATOS_FOUR_FIELD_COUPLED_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_FOUR_FIELD_COUPLED_NEWTON_RAPHSON_STRATEGY

// System includes
#include <iostream>
#include <vector>
#include <array>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos
{

/**
 * @class FourFieldCoupledNewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief Newton Raphson strategy for coupled temperature-concentration-structure-phase problems
 * @details This strategy solves the four-field coupled problem monolithically within each iteration
 * @author Your Name
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver>
class FourFieldCoupledNewtonRaphsonStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    
    KRATOS_CLASS_POINTER_DEFINITION(FourFieldCoupledNewtonRaphsonStrategy);
    
    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    
    // Field indices for organization
    enum FieldIndex {
        TEMPERATURE_FIELD = 0,
        CONCENTRATION_FIELD = 1,
        STRUCTURE_FIELD = 2,
        PHASE_FIELD = 3,
        NUM_FIELDS = 4
    };
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * @brief Constructor with all necessary components
     */
    FourFieldCoupledNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, MoveMeshFlag),
          mpScheme(pScheme),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions),
          mMaxIterationNumber(MaxIterations),
          mInitializeWasPerformed(false)
    {
        KRATOS_TRY;
        
        // Initialize builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSolver));
        
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
        
        SetEchoLevel(1);
        this->SetRebuildLevel(2);
        
        // Initialize system matrices and vectors
        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();
        
        // Initialize coupling parameters
        InitializeCouplingParameters();
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief Destructor
     */
    ~FourFieldCoupledNewtonRaphsonStrategy() override
    {
        Clear();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief Initialize the strategy
     */
    void Initialize() override
    {
        KRATOS_TRY;
        
        if (!mInitializeWasPerformed)
        {
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TConvergenceCriteriaType::Pointer p_convergence_criteria = mpConvergenceCriteria;
            
            // Initialize scheme
            if (!p_scheme->SchemeIsInitialized())
                p_scheme->Initialize(BaseType::GetModelPart());
            
            // Initialize elements for all fields
            if (!p_scheme->ElementsAreInitialized())
                p_scheme->InitializeElements(BaseType::GetModelPart());
            
            // Initialize conditions for all fields
            if (!p_scheme->ConditionsAreInitialized())
                p_scheme->InitializeConditions(BaseType::GetModelPart());
            
            // Initialize convergence criteria
            if (!p_convergence_criteria->IsInitialized())
                p_convergence_criteria->Initialize(BaseType::GetModelPart());
            
            // Initialize field-specific variables
            InitializeFieldVariables();
            
            mInitializeWasPerformed = true;
        }
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief Main solution step with coupled iteration
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();
        
        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;
        
        // Initialize iteration
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        
        // Store initial values for all fields
        StoreInitialFieldValues();
        
        // Initialize non-linear iteration
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        
        bool is_converged = false;
        
        // Main iteration loop
        while (!is_converged && iteration_number <= mMaxIterationNumber)
        {
            // Update iteration number
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            // Build coupled system matrix and RHS
            BuildCoupledSystem(rA, rb, r_model_part, p_scheme, p_builder_and_solver);
            
            // Apply coupling terms
            ApplyCouplingTerms(rA, rb, r_model_part);
            
            // Solve linear system
            TSparseSpace::SetToZero(rDx);
            p_builder_and_solver->SystemSolve(rA, rDx, rb);
            
            // Update all fields
            UpdateCoupledFields(rDx, r_model_part, p_scheme, p_builder_and_solver);
            
            // Check convergence for all fields
            is_converged = CheckCoupledConvergence(r_model_part, r_dof_set, rA, rDx, rb);
            
            // Finalize iteration
            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
            
            if (BaseType::GetEchoLevel() > 0)
            {
                PrintIterationInfo(iteration_number, is_converged);
            }
            
            iteration_number++;
        }
        
        // Check if maximum iterations exceeded
        if (iteration_number > mMaxIterationNumber && !is_converged)
        {
            KRATOS_WARNING("FourFieldCoupledStrategy") 
                << "Maximum iterations (" << mMaxIterationNumber << ") exceeded!" << std::endl;
        }
        
        // Calculate reactions if needed
        if (mCalculateReactionsFlag)
        {
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);
        }
        
        return is_converged;
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief Clear internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;
        
        if (mpBuilderAndSolver != nullptr)
        {
            mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
            mpBuilderAndSolver->Clear();
        }
        
        if (mpA != nullptr) TSparseSpace::Clear(mpA);
        if (mpDx != nullptr) TSparseSpace::Clear(mpDx);
        if (mpb != nullptr) TSparseSpace::Clear(mpb);
        
        if (mpScheme != nullptr) mpScheme->Clear();
        
        mInitializeWasPerformed = false;
        
        KRATOS_CATCH("");
    }
    
    ///@}
    ///@name Access
    ///@{
    
    typename TSchemeType::Pointer GetScheme() { return mpScheme; }
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver() { return mpBuilderAndSolver; }
    
    void SetCouplingCoefficient(FieldIndex field1, FieldIndex field2, double value)
    {
        mCouplingMatrix[field1][field2] = value;
    }
    
    double GetCouplingCoefficient(FieldIndex field1, FieldIndex field2) const
    {
        return mCouplingMatrix[field1][field2];
    }
    
    ///@}
    
protected:
    ///@name Protected member Variables
    ///@{
    
    typename TSchemeType::Pointer mpScheme;
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;
    
    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;
    
    bool mReformDofSetAtEachStep;
    bool mCalculateReactionsFlag;
    unsigned int mMaxIterationNumber;
    bool mInitializeWasPerformed;
    
    // Coupling coefficients matrix
    std::array<std::array<double, NUM_FIELDS>, NUM_FIELDS> mCouplingMatrix;
    
    // Field-specific convergence tolerances
    std::array<double, NUM_FIELDS> mFieldTolerance;
    
    // Storage for field values
    std::vector<double> mTemperatureValues;
    std::vector<double> mConcentrationValues;
    std::vector<double> mDisplacementValues;
    std::vector<double> mPhaseValues;
    
    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * @brief Initialize coupling parameters
     */
    void InitializeCouplingParameters()
    {
        // Initialize coupling matrix (symmetric)
        for (int i = 0; i < NUM_FIELDS; ++i)
        {
            for (int j = 0; j < NUM_FIELDS; ++j)
            {
                mCouplingMatrix[i][j] = 0.0;
            }
        }
        
        // Set default coupling coefficients (adjust based on your physics)
        // Temperature-Structure coupling (thermal expansion)
        mCouplingMatrix[TEMPERATURE_FIELD][STRUCTURE_FIELD] = 1.0e-5;
        mCouplingMatrix[STRUCTURE_FIELD][TEMPERATURE_FIELD] = 1.0e-5;
        
        // Temperature-Phase coupling
        mCouplingMatrix[TEMPERATURE_FIELD][PHASE_FIELD] = 1.0;
        mCouplingMatrix[PHASE_FIELD][TEMPERATURE_FIELD] = 1.0;
        
        // Concentration-Phase coupling
        mCouplingMatrix[CONCENTRATION_FIELD][PHASE_FIELD] = 0.5;
        mCouplingMatrix[PHASE_FIELD][CONCENTRATION_FIELD] = 0.5;
        
        // Structure-Phase coupling (transformation strain)
        mCouplingMatrix[STRUCTURE_FIELD][PHASE_FIELD] = 0.1;
        mCouplingMatrix[PHASE_FIELD][STRUCTURE_FIELD] = 0.1;
        
        // Set field tolerances
        mFieldTolerance[TEMPERATURE_FIELD] = 1.0e-6;
        mFieldTolerance[CONCENTRATION_FIELD] = 1.0e-6;
        mFieldTolerance[STRUCTURE_FIELD] = 1.0e-6;
        mFieldTolerance[PHASE_FIELD] = 1.0e-6;
    }
    
    /**
     * @brief Initialize field-specific variables
     */
    void InitializeFieldVariables()
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        const unsigned int num_nodes = r_model_part.NumberOfNodes();
        
        mTemperatureValues.resize(num_nodes);
        mConcentrationValues.resize(num_nodes);
        mDisplacementValues.resize(num_nodes * 3); // 3D displacement
        mPhaseValues.resize(num_nodes);
    }
    
    /**
     * @brief Store initial field values
     */
    void StoreInitialFieldValues()
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        
        unsigned int index = 0;
        for (auto& node : r_model_part.Nodes())
        {
            // Store temperature
            if (node.HasDofFor(TEMPERATURE))
                mTemperatureValues[index] = node.FastGetSolutionStepValue(TEMPERATURE);
            
            // Store concentration
            if (node.HasDofFor(CONCENTRATION))
                mConcentrationValues[index] = node.FastGetSolutionStepValue(CONCENTRATION);
            
            // Store displacement
            if (node.HasDofFor(DISPLACEMENT_X))
            {
                mDisplacementValues[index*3] = node.FastGetSolutionStepValue(DISPLACEMENT_X);
                mDisplacementValues[index*3+1] = node.FastGetSolutionStepValue(DISPLACEMENT_Y);
                mDisplacementValues[index*3+2] = node.FastGetSolutionStepValue(DISPLACEMENT_Z);
            }
            
            // Store phase field
            if (node.HasDofFor(PHASE_FIELD))
                mPhaseValues[index] = node.FastGetSolutionStepValue(PHASE_FIELD);
            
            index++;
        }
    }
    
    /**
     * @brief Build the coupled system matrix and RHS
     */
    void BuildCoupledSystem(TSystemMatrixType& rA, TSystemVectorType& rb,
                           ModelPart& rModelPart,
                           typename TSchemeType::Pointer pScheme,
                           typename TBuilderAndSolverType::Pointer pBuilderAndSolver)
    {
        KRATOS_TRY;
        
        // Clear system
        TSparseSpace::SetToZero(rA);
        TSparseSpace::SetToZero(rb);
        
        // Build the basic system
        pBuilderAndSolver->Build(pScheme, rModelPart, rA, rb);
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief Apply coupling terms to system matrix and RHS
     */
    void ApplyCouplingTerms(TSystemMatrixType& rA, TSystemVectorType& rb,
                           ModelPart& rModelPart)
    {
        KRATOS_TRY;
        
        // This is where you add the coupling terms to the system matrix
        // The implementation depends on your specific physics
        
        for (auto& element : rModelPart.Elements())
        {
            // Get element nodes
            auto& geometry = element.GetGeometry();
            const unsigned int num_nodes = geometry.size();
            
            // Local coupling matrices
            Matrix local_coupling(num_nodes * NUM_FIELDS, num_nodes * NUM_FIELDS, 0.0);
            Vector local_rhs(num_nodes * NUM_FIELDS, 0.0);
            
            // Calculate coupling contributions
            CalculateElementCouplingContribution(element, local_coupling, local_rhs);
            
            // Assemble to global system
            AssembleCouplingContribution(rA, rb, element, local_coupling, local_rhs);
        }
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief Calculate element coupling contribution
     */
    void CalculateElementCouplingContribution(Element& rElement,
                                             Matrix& rCouplingMatrix,
                                             Vector& rCouplingRHS)
    {
        // Implementation specific to your physics
        // This should calculate the coupling terms between fields
        
        auto& geometry = rElement.GetGeometry();
        const unsigned int num_nodes = geometry.size();
        
        // Example: thermal-mechanical coupling
        for (unsigned int i = 0; i < num_nodes; ++i)
        {
            for (unsigned int j = 0; j < num_nodes; ++j)
            {
                // Temperature-displacement coupling
                double coupling_coeff = mCouplingMatrix[TEMPERATURE_FIELD][STRUCTURE_FIELD];
                
                // Add coupling terms (simplified example)
                unsigned int temp_index = i * NUM_FIELDS + TEMPERATURE_FIELD;
                unsigned int disp_index = j * NUM_FIELDS + STRUCTURE_FIELD;
                
                rCouplingMatrix(temp_index, disp_index) += coupling_coeff;
                rCouplingMatrix(disp_index, temp_index) += coupling_coeff;
            }
        }
    }
    
    /**
     * @brief Assemble coupling contribution to global system
     */
    void AssembleCouplingContribution(TSystemMatrixType& rA,
                                     TSystemVectorType& rb,
                                     Element& rElement,
                                     const Matrix& rLocalMatrix,
                                     const Vector& rLocalVector)
    {
        // Get equation IDs for the element
        Element::EquationIdVectorType equation_ids;
        rElement.EquationIdVector(equation_ids, BaseType::GetModelPart().GetProcessInfo());
        
        // Assemble local contributions to global system
        for (unsigned int i = 0; i < equation_ids.size(); ++i)
        {
            unsigned int row_id = equation_ids[i];
            
            for (unsigned int j = 0; j < equation_ids.size(); ++j)
            {
                unsigned int col_id = equation_ids[j];
                rA(row_id, col_id) += rLocalMatrix(i, j);
            }
            
            rb[row_id] += rLocalVector[i];
        }
    }
    
    /**
     * @brief Update all coupled fields
     */
    void UpdateCoupledFields(TSystemVectorType& rDx,
                           ModelPart& rModelPart,
                           typename TSchemeType::Pointer pScheme,
                           typename TBuilderAndSolverType::Pointer pBuilderAndSolver)
    {
        KRATOS_TRY;
        
        // Update solution using the scheme
        pScheme->Update(rModelPart, pBuilderAndSolver->GetDofSet(), *mpA, rDx, *mpb);
        
        // Apply any field-specific constraints or corrections
        ApplyFieldConstraints(rModelPart);
        
        // Move mesh if needed (for structural field)
        if (BaseType::MoveMeshFlag())
        {
            BaseType::MoveMesh();
        }
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief Apply field-specific constraints
     */
    void ApplyFieldConstraints(ModelPart& rModelPart)
    {
        // Apply any physical constraints
        // For example: phase field should be between 0 and 1
        for (auto& node : rModelPart.Nodes())
        {
            if (node.HasDofFor(PHASE_FIELD))
            {
                double& phase = node.FastGetSolutionStepValue(PHASE_FIELD);
                phase = std::max(0.0, std::min(1.0, phase));
            }
        }
    }
    
    /**
     * @brief Check convergence for all fields
     */
    bool CheckCoupledConvergence(ModelPart& rModelPart,
                                DofsArrayType& rDofSet,
                                TSystemMatrixType& rA,
                                TSystemVectorType& rDx,
                                TSystemVectorType& rb)
    {
        // Use the convergence criteria
        bool converged = mpConvergenceCriteria->PostCriteria(rModelPart, rDofSet, rA, rDx, rb);
        
        // Additional field-specific convergence checks if needed
        if (converged)
        {
            converged = CheckFieldConvergence(rModelPart);
        }
        
        return converged;
    }
    
    /**
     * @brief Check field-specific convergence
     */
    bool CheckFieldConvergence(ModelPart& rModelPart)
    {
        double temp_norm = 0.0, conc_norm = 0.0, disp_norm = 0.0, phase_norm = 0.0;
        unsigned int index = 0;
        
        for (auto& node : rModelPart.Nodes())
        {
            // Temperature convergence
            if (node.HasDofFor(TEMPERATURE))
            {
                double temp_diff = node.FastGetSolutionStepValue(TEMPERATURE) - mTemperatureValues[index];
                temp_norm += temp_diff * temp_diff;
            }
            
            // Concentration convergence
            if (node.HasDofFor(CONCENTRATION))
            {
                double conc_diff = node.FastGetSolutionStepValue(CONCENTRATION) - mConcentrationValues[index];
                conc_norm += conc_diff * conc_diff;
            }
            
            // Displacement convergence
            if (node.HasDofFor(DISPLACEMENT_X))
            {
                double dx_diff = node.FastGetSolutionStepValue(DISPLACEMENT_X) - mDisplacementValues[index*3];
                double dy_diff = node.FastGetSolutionStepValue(DISPLACEMENT_Y) - mDisplacementValues[index*3+1];
                double dz_diff = node.FastGetSolutionStepValue(DISPLACEMENT_Z) - mDisplacementValues[index*3+2];
                disp_norm += dx_diff*dx_diff + dy_diff*dy_diff + dz_diff*dz_diff;
            }
            
            // Phase field convergence
            if (node.HasDofFor(PHASE_FIELD))
            {
                double phase_diff = node.FastGetSolutionStepValue(PHASE_FIELD) - mPhaseValues[index];
                phase_norm += phase_diff * phase_diff;
            }
            
            index++;
        }
        
        // Check tolerances
        bool temp_converged = std::sqrt(temp_norm) < mFieldTolerance[TEMPERATURE_FIELD];
        bool conc_converged = std::sqrt(conc_norm) < mFieldTolerance[CONCENTRATION_FIELD];
        bool disp_converged = std::sqrt(disp_norm) < mFieldTolerance[STRUCTURE_FIELD];
        bool phase_converged = std::sqrt(phase_norm) < mFieldTolerance[PHASE_FIELD];
        
        return temp_converged && conc_converged && disp_converged && phase_converged;
    }
    
    /**
     * @brief Print iteration information
     */
    void PrintIterationInfo(unsigned int IterationNumber, bool IsConverged)
    {
        if (IsConverged)
        {
            KRATOS_INFO("FourFieldCoupledStrategy") 
                << "Converged in " << IterationNumber << " iterations" << std::endl;
        }
        else
        {
            KRATOS_INFO("FourFieldCoupledStrategy") 
                << "Iteration " << IterationNumber << " completed" << std::endl;
        }
    }
    
    ///@}
    
}; // Class FourFieldCoupledNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_FOUR_FIELD_COUPLED_NEWTON_RAPHSON_STRATEGY