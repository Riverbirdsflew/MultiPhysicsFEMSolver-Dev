// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_EULERIAN_THERMAL_CONDUCTION_ELEMENT_INCLUDED )
#define  KRATOS_EULERIAN_THERMAL_CONDUCTION_ELEMENT_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/phase_transformation_law.h"



namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(THERMAL_CONDUCTION_APPLICATION) EulerianThermalConductionElement
    : public Element
{
public:
    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EulerianThermalConductionElement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default constructor.

    EulerianThermalConductionElement() : Element()
    {
    }

    EulerianThermalConductionElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    EulerianThermalConductionElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~EulerianThermalConductionElement() {};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& ThisNodes, 
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<EulerianThermalConductionElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }
    
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<EulerianThermalConductionElement>(NewId, pGeom, pProperties);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //返回单元类型
    std::string Info() const override
    {
        return "EulerianThermalConductionElement #";
    }
    
    //返回单元名称
    std::string GetElementTypeAsString() const override
    {
        std::string BaseElementName = Element::GetElementTypeAsString();
        
        return "EulerianThermCond" + BaseElementName;
        }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        double theta;
        double dt_inv;
        double lumping_factor; //用于平均的参数
        double conductivity;
        double specific_heat;
        double density;
        //多相材料的质量分数
        
        //当前步和上一步的温度
        array_1d<double,TNumNodes> t_now;
        array_1d<double,TNumNodes> t_old;
        
        //当前步和上一步的相分数，七种相
        array_1d< array_1d<double,7 >, TNumNodes>  mass_fraction_current;
        array_1d< array_1d<double,7 >, TNumNodes> mass_fraction_previous;

        array_1d<double,TNumNodes> volumetric_source;
        array_1d<double,TNumNodes> austenize_latent;
        array_1d<double,TNumNodes> eutectoid_latent;
        array_1d<double,TNumNodes> bainite_latent;
        array_1d<double,TNumNodes> martensite_latent;
    };

    void InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume);

    void GetNodalValues(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) const;

    /* phase transformation part */
	/* return the austenite mass fraction increment */
	virtual int GetAusteniteMassFrcInc(Kratos::PhaseTransformationLaw::Parameters &rValues, double *mfinc);
	/* return the eutectoid phase mass fraction increment */
	virtual int GetEutectoidPhaseMassFrcInc(Kratos::PhaseTransformationLaw::Parameters &rValues, double *mfinc);
	/* return the bainite mass fraction increment */
	virtual int GetBainiteMassFrcInc(Kratos::PhaseTransformationLaw::Parameters &rValues, double *mfinc);
	/* return the martensite mass fraction increment */
	virtual int GetMartensiteMassFrcInc(Kratos::PhaseTransformationLaw::Parameters &rValues, double *mfinc);
	/* update the incubation, carbon content and mass fraction of all microstructures */
	virtual int UpdateMassFraction(Kratos::PhaseTransformationLaw::Parameters &rValues, double *mfinc);

    // Member Variables


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

        // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }


}; // Class EulerianThermalConductionElement

} // namespace Kratos.

#endif // KRATOS_EULERIAN_THERMAL_CONDUCTION_ELEMENT_INCLUDED  defined
