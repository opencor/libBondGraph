/*******************************************************************************

Copyright (C) The University of Auckland

OpenCOR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenCOR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://gnu.org/licenses>.

*******************************************************************************/

#pragma once
#include "RCP.h"
#include "export.h"
#include "thirdparty/json.hpp"
#include <map>
#include <unordered_map>
#include <vector>

namespace SymEngine {
class DenseMatrix;
class Basic;
struct Comparator {
  bool operator()(const RCPLIB::RCP<const Basic> &x,
                  const RCPLIB::RCP<const Basic> &y) const;
};
} // namespace SymEngine

namespace units {
class precise_unit;
}

namespace BG {

enum EXPORTED ElementType {
  eSource,   //!< Source type
  ePassive,  //!< Passive type
  eJunction, //!< Junction type
  eBond      //!< BondInterface type
};

//! Component groups
enum EXPORTED ComponentGroup {
  eU,  //!< Source element (Se,Sf)
  eS,  //!< Storage element (I,C)
  eR,  //!< Dissipative element (R,Re)
  eJ,  //!< Junction with power conservation (J,GY,TR)
  ePH, //!< PortInterface Hamiltonians
  nogroup
};

//! Define the type of storage bond to complement the BondGroup
enum EXPORTED BondStorageType {
  eC, //!< Capacitor storage type
  eI  //!< Inductor storage type
};

// Transformater and gyrator in junction group.
//! Passive type
enum EXPORTED PassiveType {
  eResistance,      //!< A resistance component
  eCapacitance,     //!< A capacitance component
  eInductance,      //!< A inertial/inductance component
  eTransformer,     //!< A tranformator
  eGyrator,         //!< A gyrator
  ePotentialSource, //!< A Source element
  eFlowSource,      //!< A Flow element
  bReaction,        //!< A Chemical reaction
  bConcentration,   //!< A Chemical concentration
  bChemostat,       //!< A Chemical chemostat
  bFlowstat,        //!< A Chemical fluxstat
  bStoichiometry,   //!< A Chemical transformer element
  eZero,            //!< O-junction
  eOne,             //!< 1-junction
  ePortHamiltonian, //!< PortInterface hamiltonian type
  eUserDefined,     //!< User defined blackbox
  notype
};

typedef std::string PhysicalDomain;

class BGElement;
class BondInterface;

// Forward declaration
class EXPORTED PortInterface {
public:
  //! Set identification string
  virtual void setId(std::string id) = 0;
  //! Get identification string
  virtual std::string getId() = 0;
  //! Get Power direction
  virtual bool receivePower() = 0;
  //!  Connect to a component.
  virtual void connect(const RCPLIB::RCP<BGElement> &inComponent) = 0;
  //!  Connect to a bond.
  virtual void connectBond(const RCPLIB::RCP<BondInterface> &inBond) = 0;
  //! Release
  virtual void release() = 0;
  // Set dofIndex
  virtual void setDofIndex(long int ix) = 0;
  // Get dofIndex
  virtual long int dofIndex() = 0;
  //!  Get the connected component.
  virtual RCPLIB::RCP<BGElement> getComponent() const = 0;
  //!  Get the connected bond.
  virtual RCPLIB::RCP<BondInterface> getBond() const = 0;
  // Set weight
  virtual void setWeight(double wt) = 0;
  //! Get weight
  virtual const double getWeight() const = 0;
};

class EXPORTED BondInterface {
public:
  //! Get the incoming bond port.
  virtual RCPLIB::RCP<PortInterface> getFromPort() const = 0;
  //! Get the outcoming bond port.
  virtual RCPLIB::RCP<PortInterface> getToPort() const = 0;
  //! Get the other end of the bond.
  virtual RCPLIB::RCP<PortInterface>
  getOtherEnd(const RCPLIB::RCP<PortInterface> inPort) const = 0;
  //! Set the identification number of the bond.
  virtual void setId(long int inId) = 0;
  //! Get the identification number of the bond.
  virtual long int getId() = 0;
  //! Connect two ports together.
  virtual void connect(const RCPLIB::RCP<PortInterface> &inFromPort,
                       const RCPLIB::RCP<PortInterface> &inToPort) = 0;
};

/*! \brief Numeric values with units and magnitude information
 *	Structure to store state, parameter values along with the SI unit
 *information SI Unit in fundamental units is stored as an 7 dimensional array.
 *The array entries determin the power The units in order are
 *Seconds,Meters,Kilogram,Kelvin,Ampere,Mole,Candela, when all the entries are
 *0, the quantity is physically dimensionless
 */
class EXPORTED Value {
public:
  std::string prefix; // Local name, value against which local parameter search
                      // will be performed
  std::string name;   // Dof dependent name
  std::string value;
  std::string units;
  bool universalConstant;
  Value(std::string _name, double v = 0.0);
  Value(std::string _name, std::string expr = "0.0");
  Value(const Value &v);

  friend std::ostream &operator<<(std::ostream &os, const Value &val);
};

// Concrete implementation of Bondgraph component, proxy components will hold a
// pointer to an instance of this
struct EXPORTED BGElementData {
  std::string mName = "";       //!< Name of the element
  std::string mId = "";         //!< Element's ID
  std::string hamiltonian = ""; //!< Hamltonian energy description (valid if the
                                //!< element is a Port Hamiltonian)
  PassiveType mElementType = notype;        //!< Type of element
  ComponentGroup mComponentGroup = nogroup; //!< Component group association.
  PhysicalDomain mDomain = "";              //! Physical Domain
  unsigned int numStates = 0;  //! Number of states for this Bondgraph element
  bool allocatedPorts = false; //! True if ports are preallocated
  bool portModified =
      false; //! True if port assignment was modified, used by junctions
  std::vector<RCPLIB::RCP<PortInterface>> mPorts =
      {}; //!< Ports connected to the component. Size depends on the type of
          //!< component.
  std::vector<std::string> constitutiveEq = {}; //!< Constitutive equations
  std::vector<int> constitutiveEqIndex =
      {}; //!< Index of constitutive equation in connectivity matrix
  std::vector<RCPLIB::RCP<Value>> mParameter =
      {}; //! Array of state and parameter values associated with this bondgraph
          //! element
};

class EXPORTED BGElement {
public:
  virtual const std::string &getName() const = 0;
  //! Set element name
  virtual void setName(const std::string &name) = 0;
  //! Set if instance is a proxy
  virtual void setProxy(bool flag) = 0;
  //! Check if this is a proxy instance
  virtual bool isProxy() = 0;
  //! Return the element ID.
  virtual std::string getId() const = 0;
  //! Set the element ID.
  virtual void setId(std::string inId) = 0;
  //! Set the element ID.
  virtual void setId(long int inId) = 0;
  //! Set the DoF ID.
  virtual void setDof(size_t inId) = 0;
  //! Get the Dof ID.
  virtual size_t getDof() = 0;
  //! Set the parent BondGraph's name
  virtual void setParent(std::string &pid) = 0;
  //! Get the parent BondGraph's name
  virtual const std::string getParent() = 0;
  //! Set the physical domain, this will not change constitutive relations or
  //! statenames and dimensions
  virtual void setDomain(PhysicalDomain domain) = 0;
  //! Get the physical domain
  virtual PhysicalDomain getDomain() = 0;
  // Get the number of states for this bond graph element
  virtual unsigned int getNumStates() = 0;
  //! Check if ports are preallocated
  virtual bool preallocatedPorts() = 0;
  //! Set portModified
  virtual void setPortModified(bool flag = true) = 0;
  //! Return the ports associated to the component.
  virtual std::vector<RCPLIB::RCP<PortInterface>> &getPorts() = 0;
  //! Return port at pin i
  virtual RCPLIB::RCP<PortInterface> &getPorts(unsigned int i) = 0;
  //! Connect the component of a port
  virtual void connect(const RCPLIB::RCP<PortInterface> &inPort,
                       int portNum = 0) = 0;
  //! Connect the component to inPort and disconnect it from outPort
  virtual void connect(const RCPLIB::RCP<PortInterface> &inPort,
                       const RCPLIB::RCP<PortInterface> &outPort) = 0;
  //! Disconnect the component of a port
  virtual void disconnect(const RCPLIB::RCP<PortInterface> &inPort) = 0;
  //! Return the component associated type
  virtual PassiveType getType() const = 0;
  //! Get constitutive equations
  virtual std::vector<std::string> &getConstitutiveEquations() = 0;
  //! Set physical dimension units, only supported for port Hamiltonian types,
  //! Transformers and Gyrators
  virtual void setSIUnit(std::string name, std::string unit) = 0;
  //! Set parameter
  virtual RCPLIB::RCP<Value> setParameter(std::string name, std::string value,
                                          std::string unit) = 0;
  //! Set value
  virtual RCPLIB::RCP<Value> setValue(std::string name, std::string value) = 0;
  //! Get Hamiltonian
  virtual std::string getHamiltonian() = 0;
  //! Set a symbol as a Universal Constant
  virtual RCPLIB::RCP<Value>
  setUniversalConstant(std::string name, double &value, std::string unit) = 0;
  //! Get the values associated with the element
  virtual std::vector<std::tuple<std::string, RCPLIB::RCP<Value>>> values() = 0;
  //! Set the PMR annotation
  virtual void setPMRAnnotation(nlohmann::json &annotation) = 0;
  //! Get the PMR annotation
  virtual nlohmann::json &getPMRAnnotation() = 0;
};

class EXPORTED ComputeEquationResults {
public:
  // false if the logic concludes that the bg is incorrect, for instance, some
  // state and power variables are evaluating to zero, it will still return the
  // equations
  bool bondGraphValidity;
  // Name map
  std::unordered_map<std::string, std::string> nameMap;
  // system of equations corresponding to each dof - dof_symbol, eqn
  std::map<RCPLIB::RCP<const SymEngine::Basic>,
           RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
      dof;
  // system of equations corresponding to each bond - bond_variable, eqn
  std::map<RCPLIB::RCP<const SymEngine::Basic>,
           RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
      bondEquations;
  // constraints on equations corresponding to each dof - dof_symbol, eqn
  std::map<RCPLIB::RCP<const SymEngine::Basic>,
           RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
      dof_constraints;
  // list of constraints
  std::vector<RCPLIB::RCP<const SymEngine::Basic>> constraints;
  // physical dimensions, value and variable type associated with every the
  // model: the last element identifies if it is a state (s), control (c) or
  // parameter (p)
  std::unordered_map<std::string, std::tuple<std::string, std::string, char>>
      physicalDimensions;
  // Annotations associated with each model variable - the relationship
  // information is available in the json
  std::unordered_map<std::string, std::vector<nlohmann::json>> annotations;
};

class EXPORTED BondGraphInterface {
public:
  /**
   * @brief Get the unique Id object
   *
   * @return std::string
   */
  virtual std::string getId() = 0;
  /**
   * @brief ID of the bondgraph that owns this bondgraph instance, empty if its
   * independent
   *
   * @return std::string
   */
  virtual std::string owner() = 0;
  /**
   * @brief Create a clone of the input bondgraph
   * A copy of the bondgraph is made by copying all elements (proxies are also
   * concretely allocated) and bonds are restablised
   *
   * @return RCPLIB::RCP<BondGraphInterface>
   */
  virtual RCPLIB::RCP<BondGraphInterface> clone() = 0;
  // Creation related methods
  /**
   * @brief Add a bondgraph element to the graph
   *
   * @param inComponent Element to be added
   */
  virtual void addComponent(const RCPLIB::RCP<BGElement> &inComponent) = 0;
  // Creation related methods
  /**
   * @brief Add a bondgraph to the graph
   *
   * @param inComponent Element to be added
   */
  virtual void
  addBondGraph(const RCPLIB::RCP<BondGraphInterface> &inComponent) = 0;
  /**
   * @brief Connect the power ports of two bondgraph elements that exist within
   * the bondgraph
   *
   * @param inOrigin The element from which power is expected to flow out
   * @param inDest   The element to which power is expected to flow in
   * @param opin Specify the port from which power should flow out: default 0
   * @param ipin Specify the port to  which power should flow in: default 0
   * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
   */
  virtual RCPLIB::RCP<BondInterface>
  connect(const RCPLIB::RCP<BGElement> &inOrigin,
          const RCPLIB::RCP<BGElement> &inDest, int opin = 0, int ipin = 0) = 0;
  /**
   * @brief Connect the power ports of two bondgraph elements that exist within
   * the bondgraph
   *
   * @param inOrigin The element from which power is expected to flow out
   * @param opin Specify the port from which power should flow out
   * @param inDest   The element to which power is expected to flow in
   * @param ipin Specify the port to  which power should flow in
   * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
   */
  virtual RCPLIB::RCP<BondInterface>
  connect(const RCPLIB::RCP<BGElement> &inOrigin, int opin,
          const RCPLIB::RCP<BGElement> &inDest, int ipin) = 0;
  /**
   * @brief Connect the power ports of two bondgraph elements that exist within
   * the bondgraph
   *
   * @param inOrigin The element from which power is expected to flow out
   * @param opin Specify the port from which power should flow out
   * @param inDest   The element to which power is expected to flow in
   * @param ipin Specify the port to  which power should flow in: default 0
   * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
   */
  virtual RCPLIB::RCP<BondInterface>
  connect(const RCPLIB::RCP<BGElement> &inOrigin, int opin,
          const RCPLIB::RCP<BGElement> &inDest) = 0;
  /**
   * @brief Connect the power ports of two bondgraph elements that exist within
   * the bondgraph
   *
   * @param inOrigin The element from which power is expected to flow out
   * @param inDest   The element to which power is expected to flow in
   * @param ipin Specify the port to  which power should flow in
   * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
   */
  virtual RCPLIB::RCP<BondInterface>
  connect(const RCPLIB::RCP<BGElement> &inOrigin,
          const RCPLIB::RCP<BGElement> &inDest, int ipin) = 0;
  /**
   * @brief Connect the power ports of two bondgraph elements that exist within
   * the bondgraph
   *
   * @param inOrigin The element from which power is expected to flow out
   * @param opin     Specify the port from which power should flow out
   * @param inDest   The element to which power is expected to flow in
   * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
   */
  virtual RCPLIB::RCP<BondInterface>
  connectInverting(const RCPLIB::RCP<BGElement> &inOrigin, int opin,
                   const RCPLIB::RCP<BGElement> &inDest) = 0;
  /**
   * @brief Connect the power ports of two bondgraph elements that exist within
   * the bondgraph
   *
   * @param inOrigin The element from which power is expected to flow out
   * @param inDest   The element to which power is expected to flow in
   * @param ipin     Specify the port from which power should flow in
   * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
   */
  virtual RCPLIB::RCP<BondInterface>
  connectInverting(const RCPLIB::RCP<BGElement> &inOrigin,
                   const RCPLIB::RCP<BGElement> &inDest, int ipin) = 0;

  /**
   * @brief Removes a bondgraph instance with id from itself
   *
   * @param id
   */
  virtual void removeBondGraph(std::string id) = 0;

  /**
   * @brief Removes a bondgraph instance with RCP handle bg
   *
   * @param bg
   */
  virtual void removeBondGraph(const RCPLIB::RCP<BondGraphInterface> &bg) = 0;

  // Bondgraph manipulation methods
  /**
   * @brief Remove a bond
   *
   * @param inBond handle to bond
   */
  virtual void removeBond(const RCPLIB::RCP<BondInterface> &inBond) = 0;
  /**
   * @brief Remove component and associated linkages
   *
   * @param inComponent component to be removed
   */
  virtual void removeComponent(const RCPLIB::RCP<BGElement> &inComponent) = 0;
  /**
   * @brief Remove all components and their associated linkages in the list
   *
   * @param inComponents std::vector of element handles that need to be removed
   */
  virtual void
  removeComponent(std::vector<RCPLIB::RCP<BGElement>> &inComponents) = 0;
  /**
   * @brief Reconnect a component to a port
   *
   * @param inComponentPort Port to which the component needs be connected
   * @param inNewOrigin     Component that needs to be connected
   */
  virtual void
  reconnectComponent(const RCPLIB::RCP<PortInterface> &inComponentPort,
                     const RCPLIB::RCP<BGElement> &inNewOrigin) = 0;

  // State equation related methods
  /**
   * @brief Compute the state equations for the bondgraph
   * @return instance of ComputeEquationResults
   */
  virtual ComputeEquationResults computeStateEquation() = 0;

  virtual ComputeEquationResults computeStateEquationNoDim() = 0;

  /**
   * @brief Compute the Port Hamiltonian of the bondgraph
   */
  virtual nlohmann::json computePortHamiltonian() = 0;

  /**
   * @brief Get the Components object
   * Returns a tuple with readablename (if not available the id), handle to it
   * and whether it is a proxy
   * @return std::vector<std::tuple<str::string,RCPLIB::RCP<BGElement> ,bool>>&
   */
  virtual std::vector<
      std::tuple<std::string, const RCPLIB::RCP<BGElement> &, bool>>
  getComponents() = 0;

  virtual const RCPLIB::RCP<BGElement> &getElementByName(std::string name) = 0;

  /**
   * @brief Get the list of source elements associated with this bondgraph,
   * entries are valid only after computeStateEquation is called
   *
   * @return std::vector<RCPLIB::RCP<BGElement> >& List of elements
   */
  virtual std::vector<RCPLIB::RCP<BGElement>> &getSources() = 0;
  /**
   * @brief Get the list of Bonds associated with this bondgraph, entries are
   * valid only after computeStateEquation is called
   *
   * @return const std::vector<RCPLIB::RCP<BondInterface>>& List of bonds
   */
  virtual const std::vector<RCPLIB::RCP<BondInterface>> &getBonds() const = 0;

  // IO related methods
  //! Serialize the bond graph in JSON form.
  virtual std::string serialize(bool inIndent = false) const = 0;
  //! Serialize the bond graph in flattend one component CellML form.
  virtual std::string serializeAsCellML(std::string modelName) = 0;
  virtual void read(std::string filename) = 0;
};

class EXPORTED ComponentRegistryInterface {
public:
  virtual void ownComponent(std::string parent, RCPLIB::RCP<BGElement> comp,
                            std::string readableName = "") = 0;
  virtual void ownBondgraph(std::string parent,
                            RCPLIB::RCP<BondGraphInterface> comp,
                            std::string readableName = "") = 0;
  // These functions do not take ownership of component memory management
  virtual void addComponent(std::string parent,
                            const RCPLIB::RCP<BGElement> &comp,
                            std::string readableName = "") = 0;
  virtual void addBondgraph(std::string parent,
                            const RCPLIB::RCP<BondGraphInterface> &comp,
                            std::string readableName = "") = 0;
  virtual void setName(const RCPLIB::RCP<BGElement> &comp,
                       std::string readableName) = 0;
  virtual void setName(const RCPLIB::RCP<BondGraphInterface> &comp,
                       std::string readableName) = 0;
  virtual void setName(std::string id, std::string readableName) = 0;
  virtual void removeComponentByName(std::string readableName) = 0;
  virtual void removeComponent(std::string id) = 0;
  virtual void removeComponentIgnoreIfNotAvailable(std::string id) = 0;
  virtual std::tuple<int, const RCPLIB::RCP<BGElement> &,
                     const RCPLIB::RCP<BondGraphInterface> &>
  getComponentByName(std::string readableName) = 0;
  virtual std::tuple<int, const RCPLIB::RCP<BGElement> &,
                     const RCPLIB::RCP<BondGraphInterface> &>
  getComponent(std::string id) = 0;
  virtual RCPLIB::RCP<nlohmann::json>
  getPropertiesByName(std::string readableName) = 0;
  virtual RCPLIB::RCP<nlohmann::json> getProperties(std::string id) = 0;
  virtual std::unordered_map<std::string, std::string> &getNameMap() = 0;
};

/**
 * @brief Get pointer to default component registry
 *
 * @return ComponentRegistryInterface&
 */
EXPORTED ComponentRegistryInterface &defaultRegistry();

/**
 * @brief Create a Registry object
 *
 */
EXPORTED void createRegistry();

/**
 * @brief Clears existing workspace, frees memory and starts afresh
 * User is responsible for saving any Changes in the objects prior to calling
 * this function
 */
EXPORTED void newWorkSpace();

/**
 * @brief Create a Bond Graph instance
 *
 * @return Reference counted pointer to the bondgraph
 */
EXPORTED RCPLIB::RCP<BondGraphInterface> createBondGraph();

#include "factorymethods.h"

/**
 * @brief Create a Effort Source instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createPotentialSource(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Flow Source instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createFlowSource(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Concentration instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createConcentration(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Transformer instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createTransformer(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Gyrator instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createGyrator(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Reaction instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createReaction(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Chemostat instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createChemostat(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Flowstat instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createFlowstat(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a User defined instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createUserDefined(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Port Hamiltonian object
 *
 * @param expr - Hamiltonian expression
 * @param states - std::vector<string> of symbols that form the states
 * @param proxy  - create a proxy of this ph instance
 * @return EXPORTED
 */
EXPORTED RCPLIB::RCP<BGElement>
createPortHamiltonian(std::string expr = "",
                      std::vector<std::string> states = {},
                      const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a One Junction instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createOneJunction(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create a Zero Junction instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement>
createZeroJunction(const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);

/**
 * @brief Create proxy instance
 *
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BGElement> createProxy(const RCPLIB::RCP<BGElement> &obj);

/**
 * @brief Create a Clone object
 *
 * @param Reference counted pointer to the element
 * @return RCPLIB::RCP<BGElement>
 */
EXPORTED RCPLIB::RCP<BGElement> createClone(const RCPLIB::RCP<BGElement> &obj);

/**
 * @brief Create a Bondgraph Element object
 *
 * @param methodName - factory method name
 * @param expr       - porthamiltonian factory method requires an expression,
 * and states other do not. For porthamiltonians's this must be a json object
 * string like {"Hamiltonian":"x^2+y^3-3.0*a","states":["x","y"]}
 * @return EXPORTED
 */
EXPORTED RCPLIB::RCP<BGElement> createBondgraphElement(std::string methodName,
                                                       std::string expr = "");
/**
 * @brief Clears all bondgraph(s), annotation(s) and elements from the registry
 * and frees associated memory User is expected to save any changes to these
 * instances prior to this call
 */
EXPORTED void newWorkSpace();

/**
 * @brief Check if the unit specification (unit) matches the base unit type
 * (target)
 * @param unit the unit spec that needs to be checked
 * @param target the physical unit against which the check needs to be carried
 * out
 * @return std::string empty if the base units match, else base units of `unit`
 * will be returned
 */

EXPORTED std::string checkUnits(std::string unit, std::string target);

/**
 * @brief Get the CellML object
 *
 * @param modelName
 * @param host
 * @param bgequations
 * @return std::map<std::string, std::string>, a map of filename , contents
 * which can be serialised by the calling logic
 */

EXPORTED std::map<std::string, std::string>
getCellML(std::string modelName, const RCPLIB::RCP<BondGraphInterface> &host,
          ComputeEquationResults &bgequations);

/**
 * @brief get the std string for a symengine basic expression
 *
 * @param bs
 * @return std::string
 */
EXPORTED std::string
symEngineExpressionToString(const RCPLIB::RCP<const SymEngine::Basic> &bs);

/**
 * @brief Get the Supported Physical Domains And Factory Methods as json object
 *
 * @return nlohmann::json
 */
EXPORTED nlohmann::json getSupportedPhysicalDomainsAndFactoryMethods();

/**
 * @brief Generate a bondgraph from json description of element composition
 *
 * @param string - stringified json
 * @return Managed pointer to BondGraphInterface
 */
EXPORTED RCPLIB::RCP<BondGraphInterface> generateBondGraph(std::string bgJson);

/**
 * @brief Generate port hamiltonian for a bondgraph described in the input json
 * string
 *
 * @param string - stringified json
 * @return nlohmann::json of the attempt
 */
EXPORTED nlohmann::json generatePortHamiltonian(std::string bgJson);

/**
 * @brief Support functions for json translation
 *
 */
EXPORTED RCPLIB::RCP<BGElement>
loadJson(const nlohmann::json &j,
         const RCPLIB::RCP<BGElementData> &proxy = RCPLIB::null);
EXPORTED void to_json(nlohmann::json &j, const RCPLIB::RCP<BGElementData> &p);
EXPORTED void to_json(nlohmann::json &j, const RCPLIB::RCP<PortInterface> &p_);
EXPORTED void to_json(nlohmann::json &j, const RCPLIB::RCP<Value> &p);
EXPORTED void to_json(nlohmann::json &j, const RCPLIB::RCP<BondInterface> &b);
EXPORTED void to_json(nlohmann::json &j, const RCPLIB::RCP<BGElement> &p_);
EXPORTED void to_json(nlohmann::json &j,
                      const RCPLIB::RCP<BondGraphInterface> &p_);
EXPORTED void from_json(const nlohmann::json &j,
                        RCPLIB::RCP<BondGraphInterface> &p_);
EXPORTED void from_json(const nlohmann::json &j, RCPLIB::RCP<Value> &p);
EXPORTED void from_json(const nlohmann::json &j, units::precise_unit &p);
EXPORTED void to_json(nlohmann::json &j, const units::precise_unit &p);

} // namespace BG