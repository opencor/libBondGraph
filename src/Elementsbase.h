#ifndef ELEMENTSBASE_H_
#define ELEMENTSBASE_H_
#include "RCP.h"
#include "friends.h"
#include "export.h"
#include <string>
#include <units.hpp>
#include <unordered_map>
#include <vector>

namespace BG {

enum ElementType
{
    eSource, //!< Source type
    ePassive, //!< Passive type
    eJunction, //!< Junction type
    eBond //!< Bond type
};

//! Component groups
enum ComponentGroup
{
    eU, //!< Source element (Se,Sf)
    eS, //!< Storage element (I,C)
    eR, //!< Dissipative element (R,Re)
    eJ, //!< Junction with power conservation (J,GY,TR)
    ePH, //!< Port Hamiltonians
    nogroup
};

//! Define the type of storage bond to complement the BondGroup
enum BondStorageType
{
    eC, //!< Capacitor storage type
    eI //!< Inductor storage type
};

// Transformater and gyrator in junction group.
//!Passive type
enum PassiveType
{
    eResistance, //!< A resistance component
    eCapacitance, //!< A capacitance component
    eInductance, //!< A inertial/inductance component
    eTransformer, //!< A tranformator
    eGyrator, //!< A gyrator
    ePotentialSource, //!< A Potential Source element
    eFlowSource, //!< A Flow element
    bReaction, //!< A Chemical reaction
    bConcentration, //!< A Chemical concentration
    eZero, //!< O-junction
    eOne, //!< 1-junction
    ePortHamiltonian, //!< Port hamiltonian type
    eUserDefined, //!< User defined blackbox
    notype
};


typedef std::string PhysicalDomain;


//Forward declaration
class Port;
class Bond;
class ElementImpl;
class Value;
class BondGraph;

/*! \brief Bond Graph BondGraphElement
    *  Encapsulates a bond graph element. 
    */

class BondGraphElementBase
{
protected:
    RCPLIB::RCP<ElementImpl> _data;
    std::string parent = "";
    bool proxy = false; //! True if this is a proxy, ElementImpl data is not owned by this instance, the owner has the same id

    std::string &mName; //!< Name of the element
    std::string &hamiltonian; //!< Hamltonian energy description (valid if the element is a Port Hamiltonian)
    std::string &mId; //!< Element's ID
    PassiveType &mElementType; //!< Type of element
    ComponentGroup &mComponentGroup; //!< Component group association.
    PhysicalDomain &mDomain; //! Physical Domain
    unsigned int &numStates; //Number of states for this Bondgraph element
    bool &allocatedPorts; //! True if ports are preallocated
    bool &portModified; //! True if port assignment was modified, used by junctions
    std::vector<RCPLIB::RCP<Port>> &mPorts; //!< Ports connected to the component. Size depends on the type of component.
    std::vector<std::string> &constitutiveEq; //!< Constitutive equations
    std::vector<int> &constitutiveEqIndex; //!< Index of constitutive equation in connectivity matrix
    std::vector<RCPLIB::RCP<Value>> &mParameter; //!Array of state and parameter values associated with this bondgraph element
    RCPLIB::RCP<BondGraphElementBase> rcpPtr; //RCP pointer for this
    size_t dofID = 0; //! Dof id  for the element
    //Map of physical domains
    DEFINE_PHYSICAL_DOMAINS
    BondGraphElementBase(RCPLIB::RCP<ElementImpl> data, std::string inId = "-1");

    //! Proxy constructor, the data is shared and not created
    BondGraphElementBase(const RCPLIB::RCP<BondGraphElementBase> &data);

public:
    virtual ~BondGraphElementBase();
    //! Get element name.
    const std::string &getName() const;
    //! Set element name
    void setName(const std::string &name);
    //! Set if instance is a proxy
    void setProxy(bool flag);
    //! Check if this is a proxy instance
    bool isProxy();
    //! Return the element ID.
    std::string getId() const;
    //! Set the element ID.
    virtual void setId(std::string inId);
    //! Set the element ID.
    virtual void setId(long int inId);
    //! Set the DoF ID.
    virtual void setDof(size_t inId);
    //! Get the Dof ID.
    virtual size_t getDof();
    //! Set the parent BondGraph's name
    virtual void setParent(std::string &pid);
    //! Get the parent BondGraph's name
    virtual const std::string getParent();
    //! Set the physical domain, this will not change constitutive relations or statenames and dimensions
    void setDomain(PhysicalDomain domain);
    //! Get the physical domain
    PhysicalDomain getDomain();
    //Get the number of states for this bond graph element
    unsigned int getNumStates();
    //! Check if ports are preallocated
    bool preallocatedPorts();
    //! Set portModified
    void setPortModified(bool flag = true);
    //! Return the ports associated to the component.
    std::vector<RCPLIB::RCP<Port>> &getPorts();
    //! Return port at pin i
    RCPLIB::RCP<Port> &getPorts(unsigned int i);
    //! Connect the component of a port
    virtual void connect(const RCPLIB::RCP<Port> &inPort);
    //! Connect the compomemy to inPort and disconnect it from outPort
    void connect(const RCPLIB::RCP<Port> &inPort, const RCPLIB::RCP<Port> &outPort);
    //! Disconnect the component of a port
    virtual void disconnect(const RCPLIB::RCP<Port> &inPort);
    //! Return the component associated type
    PassiveType getType() const;
    //! Return the component associated group.
    ComponentGroup getComponentGroup() const;
    //! Set constitutive equations index
    virtual void setConstitutiveEqIndex(int index, int eqNo = 0);
    //! Get constitutive equations index
    int getConstitutiveEqIndex(int eqNo = 0);
    //! Get constitutive equations
    virtual std::vector<std::string> &getConstitutiveEquations();
    //! Set physical dimension units, only supported for port Hamiltonian types, Transformers and Gyrators
    virtual void setSIUnit(std::string name, std::string unit);
    //! Set parameter
    template<typename T>
    RCPLIB::RCP<Value> setParameter(std::string name, T value, std::string unit);
    //! Set value
    template<typename T>
    RCPLIB::RCP<Value> setValue(std::string name, T value);
    //! Set parameter
    virtual RCPLIB::RCP<Value> setParameter(std::string name, std::string value, std::string unit);
    //! Set value
    virtual RCPLIB::RCP<Value> setValue(std::string name, std::string value);
    //! Get Hamiltonian
    virtual std::string getHamiltonian();
    //! Set a symbol as a Universal Constant
    virtual RCPLIB::RCP<Value> setUniversalConstant(std::string name, double &value, std::string unit);
    //! Get the values associated with the element
    virtual std::vector<std::tuple<std::string, RCPLIB::RCP<Value>>> values();

    friend RCPLIB::RCP<BondGraphElementBase> createProxy(const RCPLIB::RCP<BondGraphElementBase> &obj);
    friend RCPLIB::RCP<BondGraphElementBase> createClone(const RCPLIB::RCP<BondGraphElementBase> &obj);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraphElementBase> &p);
    friend void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraph> &p);
};

std::string replaceAll(std::string str, const std::string &from, const std::string &to);
} // namespace BG

#endif