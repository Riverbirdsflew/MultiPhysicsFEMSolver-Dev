//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

// System includes


// External includes


// Include Base h
#include "concentration_face_condition.h"
#include "utilities/integration_utilities.h"
#include "includes/table.h"

namespace Kratos
{

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


/**
 * Constructor using Geometry
 */
ConcentrationFace::ConcentrationFace(IndexType NewId, Geometry< Node >::Pointer pGeometry):
    Condition(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
ConcentrationFace::ConcentrationFace(
	IndexType NewId, 
	Geometry< Node >::Pointer pGeometry, 
	Properties::Pointer pProperties):
    Condition(NewId, pGeometry, pProperties) 
{
}


/**
 * Destructor
 */
ConcentrationFace::~ConcentrationFace() { }

///@}
///@name Operators
///@{


///@}
///@name Operations
///@{

/**
 * CONDITIONS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

// Public Operations //////////////////////////////////////////////////////////

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param ThisNodes: the nodes of the new condition
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer ConcentrationFace::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcentrationFace>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer ConcentrationFace::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConcentrationFace>(NewId, pGeom, pProperties);
}



/**
 * this determines the condition equation ID vector for all conditional
 * DOFs
 * @param rResult: the condition equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ConcentrationFace::EquationIdVector(
	EquationIdVectorType& rResult, 
	const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the equation ids. vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rResult.size() != n_nodes) {
        rResult.resize(n_nodes, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    for (unsigned int i = 0; i < n_nodes; ++i){
        rResult[i] = r_geometry[i].GetDof(CONCENTRATION).EquationId();
    }

    KRATOS_CATCH("")
}

/**
 * determines the condition equation list of DOFs
 * @param ConditionDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void ConcentrationFace::GetDofList(
	DofsVectorType& rConditionalDofList, 
	const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the DOFs vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rConditionalDofList.size() != n_nodes){
        rConditionalDofList.resize(n_nodes);
    }

    // Fill the DOFs vector from the condition nodes
    for (unsigned int i = 0; i < n_nodes; ++i){
        rConditionalDofList[i] = r_geometry[i].pGetDof(CONCENTRATION);
    }

    KRATOS_CATCH("")
}

/**
 * CONDITIONS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
 * they can be managed internally with a private method to do the same calculations
 * only once: MANDATORY
 */

void ConcentrationFace::SetIntegrationWeight(
    const IndexType IntegrationPointIndex,
    const typename GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const Vector &rJacobianDeterminantsVector,
    ConditionDataStruct &rData)
{
    rData.Weight = rJacobianDeterminantsVector[IntegrationPointIndex] * rIntegrationPoints[IntegrationPointIndex].Weight();
}

/**
 * this is called during the assembling process in order
 * to calculate all condition contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void ConcentrationFace::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix only
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ConcentrationFace::CalculateLeftHandSide(
	MatrixType& rLeftHandSideMatrix, 
	const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

    // Check (and resize) LHS matrix
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_nodes,n_nodes);

    // Declare Gauss pt. data container
    ConditionDataStruct gauss_pt_data;
    this->FillConditionDataStructure(rCurrentProcessInfo, gauss_pt_data);

    // Get geometry data
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    Vector gauss_pts_J_det = ZeroVector(n_gauss);
    r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Gauss pts. loop
    for (unsigned int g = 0; g < n_gauss; g++) {
        gauss_pt_data.N = row(N_container, g);
        SetIntegrationWeight(g, r_gauss_pts, gauss_pts_J_det, gauss_pt_data);
        this->AddIntegrationPointLHSContribution(rLeftHandSideMatrix, gauss_pt_data);
    }

    KRATOS_CATCH("")
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector only
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ConcentrationFace::CalculateRightHandSide(
	VectorType& rRightHandSideVector, 
	const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

    // Check (and resize) RHS vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes,false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(n_nodes);

    // Declare Gauss pt. data container
    ConditionDataStruct gauss_pt_data;
    this->FillConditionDataStructure(rCurrentProcessInfo, gauss_pt_data);

    // Get geometry data
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    Vector gauss_pts_J_det = ZeroVector(n_gauss);
    r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Gauss pts. loop
    for (unsigned int g = 0; g < n_gauss; g++) {
        gauss_pt_data.N = row(N_container, g);
        SetIntegrationWeight(g, r_gauss_pts, gauss_pts_J_det, gauss_pt_data);
        this->AddIntegrationPointRHSContribution(rRightHandSideVector, gauss_pt_data);
    }

    KRATOS_CATCH("")
}

inline GeometryData::IntegrationMethod ConcentrationFace::GetIntegrationMethod() const
{
    return IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(this->GetGeometry());
}

void ConcentrationFace::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> > &rVariable,
    std::vector<array_1d<double,3> > &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    if (rVariable == NORMAL) {
        const auto &r_geometry = this->GetGeometry();
        const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
        for (unsigned int g = 0; g < n_gauss; ++g) {
            rValues[g] = r_geometry.Normal(r_gauss_pts[g].Coordinates());
        }
    } else {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const ConcentrationFace* const_this = static_cast< const ConcentrationFace* >(this);
        rValues[0] = const_this->GetValue(rVariable);
        // Copy the values to the different gauss points
        for (unsigned int g = 1; g < n_gauss; ++g) {
            noalias(rValues[g]) = rValues[0];
        }
    }
}

void ConcentrationFace::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ConcentrationFace* const_this = static_cast< const ConcentrationFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        rValues[g] = rValues[0];
    }
}

void ConcentrationFace::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ConcentrationFace* const_this = static_cast< const ConcentrationFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        noalias(rValues[g]) = rValues[0];
    }
}

void ConcentrationFace::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ConcentrationFace* const_this = static_cast< const ConcentrationFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        noalias(rValues[g]) = rValues[0];
    }
}

void ConcentrationFace::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ConcentrationFace* const_this = static_cast< const ConcentrationFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        noalias(rValues[g]) = rValues[0];
    }
}

///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a string.

std::string ConcentrationFace::Info() const {
    std::stringstream buffer;
    buffer << "ConcentrationFace #" << Id();
    return buffer.str();
}

/// Print information about this object.

void ConcentrationFace::PrintInfo(std::ostream& rOStream) const {
    rOStream << "ConcentrationFace #" << Id();
}

/// Print object's data.

void ConcentrationFace::PrintData(std::ostream& rOStream) const {
	rOStream << "ConcentrationFace #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Finite element functions ///////////////////////////////////////////////////

void ConcentrationFace::AddIntegrationPointRHSContribution(
    VectorType& rRightHandSideVector,
    const ConditionDataStruct &rData)
{
    const double gauss_pt_unknown = rData.GaussPointUnknown();
    const double gauss_pt_flux = rData.GaussPointFaceConcentrationFlux();
    const double aux_concen_rhs = rData.ConcenConvectionCoefficient * (gauss_pt_unknown - rData.AmbientConcentration);
    for (unsigned int i = 0; i < (this->GetGeometry()).PointsNumber(); ++i) {
        // Add external face concentration flux contribution
        rRightHandSideVector[i] += rData.N(i) * gauss_pt_flux * rData.Weight;

        // Add ambient concentration convection contribution
        rRightHandSideVector[i] -= rData.N(i) * aux_concen_rhs * rData.Weight;
    }
}

void ConcentrationFace::AddIntegrationPointLHSContribution(
    MatrixType& rLeftHandSideMatrix,
    const ConditionDataStruct &rData)
{
    const unsigned int n_nodes = (this->GetGeometry()).PointsNumber();

    for (unsigned int i = 0; i < n_nodes; ++i) {
        for (unsigned int j = 0; j < n_nodes; ++j) {

            // Ambient temperature convection contribution
            rLeftHandSideMatrix(i,j) += rData.N(i) * rData.ConcenConvectionCoefficient * rData.N(j) * rData.Weight;
        }
    }
}

void ConcentrationFace::FillConditionDataStructure(
    const ProcessInfo &rCurrentProcessInfo,
    ConditionDataStruct &rData)
{
    // Set user-defined variables values in data container
    const auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    rData.UnknownValues.resize(n_nodes, false);
    rData.FaceConcentrationFluxValues.resize(n_nodes, false);

    for (unsigned int i = 0; i < n_nodes; ++i) {
        rData.UnknownValues[i] = r_geometry[i].FastGetSolutionStepValue(CONCENTRATION);
        rData.FaceConcentrationFluxValues[i] = r_geometry[i].FastGetSolutionStepValue(FACE_CONCENTRATION_FLUX);
    }

    // Check (and resize) and fill data container arrays
    if (rData.UnknownValues.size() != n_nodes) {
        rData.UnknownValues.resize(n_nodes, false);
    }
    if (rData.FaceConcentrationFluxValues.size() != n_nodes) {
        rData.FaceConcentrationFluxValues.resize(n_nodes, false);
    }

    // Fill data container values from properties
    // Const reference is required to have thread-safe access
    const auto &r_prop = this->GetProperties();
    //rData.Emissivity = r_prop[EMISSIVITY];
    rData.AmbientConcentration = r_prop[AMBIENT_CONCENTRATION];
    rData.ConcenConvectionCoefficient = r_prop[CONCEN_CONVECTION_COEFFICIENT];
}

///@}
///@name Friends
///@{

///@}

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

///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

///@}
///@name Serialization
///@{

ConcentrationFace::ConcentrationFace():
    Condition()
{
}

void ConcentrationFace::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );

    // List
    // To be completed with the class member list
}

void ConcentrationFace::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );

    // List
    // To be completed with the class member list
}

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
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


} // namespace Kratos.
