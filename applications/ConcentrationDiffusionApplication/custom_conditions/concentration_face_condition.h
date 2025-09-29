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

#if !defined(KRATOS_CONCENTRATION_FACE_H_INCLUDED )
#define  KRATOS_CONCENTRATION_FACE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

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

//ConcentrationFace类继承自Condition类，宏定义KRATOS_API控制类的导出方式
class KRATOS_API(CONCENTRATION_DIFFUSION_APPLICATION) ConcentrationFace: public Condition
{
public:

    

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ConcentrationFace
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConcentrationFace);
	
	/**
     * @brief Gauss pt. data structure
     * Auxiliar data structure to pass the Gauss pt. data
     */
    struct ConditionDataStruct
    {
        double Weight;              // Gauss point weight
        array_1d<double, 3> Normal; // Condition normal
        Vector N;                   // Gauss point shape functions values

        double AmbientConcentration;    // Ambient concentration value
        double ConcenConvectionCoefficient; // Ambient concentration convection coefficient
        Vector UnknownValues;         // Previous iteration unknown values
        Vector FaceConcentrationFluxValues;    // Nodal face matter concentration flux values

        double inline GaussPointUnknown() const
        {
            return InterpolateInGaussPoint(UnknownValues);
        }

        double inline GaussPointFaceConcentrationFlux() const
        {
            return InterpolateInGaussPoint(FaceConcentrationFluxValues);
        }

        double inline InterpolateInGaussPoint(const Vector &rNodalValues) const
        {
            double gauss_pt_val = 0.0;
            for (unsigned int i = 0; i < N.size(); ++i) {
                gauss_pt_val += N[i] * rNodalValues[i];
            }
            return gauss_pt_val;
        }
    };

    typedef Condition::MatrixType MatrixType;
    typedef Condition::VectorType VectorType;

    //~ /// Stefan Boltzmann constant for radiation in SI units: [W / (m^2 K^4)].
    //~ constexpr static double StefanBoltzmann = 5.67e-8;//热辐射，扩散不需要
    
    
    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    ConcentrationFace(
		IndexType NewId,
		Geometry< Node >::Pointer pGeometry);
		
    ConcentrationFace(
		IndexType NewId,
		Geometry< Node >::Pointer pGeometry,
		Properties::Pointer pProperties);		

    /**
     * Destructor
     */
    ~ConcentrationFace() override;

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

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
		IndexType NewId, 
		NodesArrayType const& ThisNodes, 
		PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
    IndexType NewId, 
    GeometryType::Pointer pGeom, 
    PropertiesType::Pointer pProperties) const override;

    /**
     * this determines the condition equation ID vector for all condition
     * DOFs
     * @param rResult: the condition equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(
		EquationIdVectorType& rResult, 
		const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * determines the condition list of DOFs
     * @param rConditionDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(
		DofsVectorType& rConditionDofList, 
		const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * CONDITIONS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all condition contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rRightHandSideVector: the condition right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(
		MatrixType& rLeftHandSideMatrix, 
		const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector: the condition right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
     
    void CalculateRightHandSide(
		VectorType& rRightHandSideVector, 
		const ProcessInfo& rCurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;


    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6 > >& rVariable,
        std::vector<array_1d<double, 6 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected LifeCycle
    ///@{
    
    // Internal default constructor for serialization
    ConcentrationFace();
    
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

    /**
     * @brief Calculates and sets the integration weight//计算积分权重
     * This function calculates the integration point weight and saves it in the provided condition data container
     * @param IntegrationPointIndex Index of current integration point
     * @param rIntegrationPoints Vector containing the integration points
     * @param rJacobianDeterminantsVector Vector containing the determinants of the Jacobian
     * @param rConditionData Condition data structure storing the integration weight
     */
    virtual void SetIntegrationWeight(
        const IndexType IntegrationPointIndex,
        const typename GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
        const Vector& rJacobianDeterminantsVector,
        ConditionDataStruct& rData);

    void AddIntegrationPointRHSContribution(
        VectorType& rRightHandSideVector,
        const ConditionDataStruct &rData);

    void AddIntegrationPointLHSContribution(
        MatrixType& rLeftHandSideMatrix,
        const ConditionDataStruct &rData);

    void FillConditionDataStructure(
        const ProcessInfo &rCurrentProcessInfo,
        ConditionDataStruct &rData);

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}

private:

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

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ConcentrationFace& operator=(ConcentrationFace const& rOther);

    /// Copy constructor.
    ConcentrationFace(ConcentrationFace const& rOther);

    ///@}

}; // Class ConcentrationFace

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_CONCENTRATION_FACE_H_INCLUDED  defined
