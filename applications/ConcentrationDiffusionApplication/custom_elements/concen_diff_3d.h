// KRATOS ___ ___  _  _  ___ ___ _  _     ___ ___ ___ ___ 
//       / __/ _ \| \| |/ __| __| \| |___|   \_ _| __| __|
//      | (_| (_) | .` | (__| _|| .` |___| |) | || _|| _| 
//       \___\___/|_|\_|\___|___|_|\_|   |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_CONCENDIFF_ELEM_3D_H_INCLUDED )
#define  KRATOS_CONCENDIFF_ELEM_3D_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "includes/serializer.h"

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

/// Short class definition.
/** A stabilized element for solving the concentration diffusion problem in 3D. @see ConcenDiff2D
*/
class ConcenDiff3D
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ConcenDiff3D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConcenDiff3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConcenDiff3D(IndexType NewId, GeometryType::Pointer pGeometry);
    ConcenDiff3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~ConcenDiff3D();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const override;

    ///Evaluates  \f$ L h s = \frac{\rho C}{\Delta t} (W, N) + (W, v. \nabla N) + \kappa (\nabla W, \nabla N)  \f$ and \f$R h s = \rho (W, Q) + \frac{\rho C}{\Delta t} (W, T^n) - L h s \ast T \f$
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
    ///Provides the global indices for each one of this element's local rows. @see NoNewtonianASGS2D
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    ///Returns a list of the element's Dofs. @see NoNewtonianASGS2D
    void GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const override;
    /// Calculates the temperature convective projection
    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

/*    double ComputeSmagorinskyViscosity(const BoundedMatrix<double, 4, 3 > & DN_DX,const double& h,const double& C, const double nu);*/



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    ConcenDiff3D() : Element()
    {
    }

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }




    ///@}
    ///@name Private Operators
    ///@{

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

}; // Class ConcenDiff3D

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_CONCENDIFF_ELEM_3D_H_INCLUDED  defined


