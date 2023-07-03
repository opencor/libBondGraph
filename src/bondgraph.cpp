#include "bondgraph.h"
#include "Bond.h"
#include "Elements.h"
#include "Elementsbase.h"
#include "Exceptions.h"
#include "Port.h"
#include "bondgraph.hpp"
#include "componentregistry.h"
#include "logging.h"
#include "manipulations.h"
#include "thirdparty/ngraph.hpp"
#include <algorithm>
#include <random>
#include <sstream>
#include <string>
#include <symengine/matrix.h>
#include <symengine/parser.h>
#include <symengine/parser/parser.h>
#include <symengine/simplify.h>
#include <symengine/solve.h>
#include <symengine/visitor.h>
#include <tuple>
#include <units.hpp>
#include <unordered_map>

#include "simplifyexplog.h"

namespace NGraph {
// Combine graphs
sGraph &operator+=(sGraph &A, const sGraph &B) {
  std::vector<sGraph::edge> b = B.edge_list();

  for (std::vector<sGraph::edge>::const_iterator t = b.begin(); t < b.end();
       t++) {
    A.insert_edge(*t);
  };

  return A;
}
// Remove graph
/**
    remove edges
*/
sGraph operator-(const sGraph &A, const sGraph &B) {
  sGraph U(A);
  std::vector<sGraph::edge> b = B.edge_list();

  for (std::vector<sGraph::edge>::const_iterator t = b.begin(); t < b.end();
       t++) {
    if (U.includes_edge(*t))
      U.remove_edge(*t);
  };

  return U;
}

} // namespace NGraph

namespace SymEngine {
// Round about to not expose symengine and its dependencies (boost etc) for
// downstream applications
bool Comparator::operator()(const RCPLIB::RCP<const Basic> &x,
                            const RCPLIB::RCP<const Basic> &y) const {
  hash_t xh = x->hash(), yh = y->hash();
  if (xh != yh)
    return xh < yh;
  if (eq(*x, *y))
    return false;
  return x->__cmp__(*y) == -1;
}
} // namespace SymEngine

namespace BG {

Value::Value(std::string _name, double v)
    : prefix(_name), name(_name), universalConstant(false) {
  value = std::to_string(v);
  units = units::to_string(units::unit_from_string("dimless"));
}

Value::Value(std::string _name, std::string expr)
    : prefix(_name), name(_name), value(expr), universalConstant(false) {
  units = units::to_string(units::unit_from_string("dimless"));
}

Value::Value(const Value &v) {
  prefix = v.prefix;
  name = v.name;
  value = v.value;
  units = v.units;
  universalConstant = v.universalConstant;
}

BondGraph::BondGraph() {
  std::string str(
      "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

  std::random_device rd;
  std::mt19937 generator(rd());

  std::shuffle(str.begin(), str.end(), generator);
  id = "ABI" + str.substr(32);
  computedCoordinates = RCPLIB::rcp(new SymEngine::DenseMatrix(0, 0));
  graph = RCPLIB::rcp(new NGraph::sGraph());
  logInfo("Bondgraph ", id, " created");
}

BondGraph::~BondGraph() {
  for (auto &c : mComponents) {
    if (!c->isProxy()) { // Let the parent do the deregistration
      // Call the Ignore version as some times spdlog closes it streams by the
      // time this destructor is called. If the component is not found, a
      // exception is thrown which will call spdlog (if debugging is enabled)
      // and cause a crash
      ComponentRegistry::defaultRegistry().removeComponentIgnoreIfNotAvailable(
          c->getId());
    }
  }
  /*
  //Remove parts
  //The child bg's will call this to remove themselves
  //Following code leads to memory leak as Teuchos crashes when memory free is
  called twice if (ownerID != "" && ownerID != id) { const auto &entry =
  ComponentRegistry::defaultRegistry().getComponent(ownerID); const auto
  &parentBg = std::get<2>(entry); if (!parentBg.is_null()) { try {
              parentBg->removeBondGraph(id);
          } catch (std::exception &ex) {
              logWarn("Exception when degregistering bondgraph with id ", id, "
  from its ownwer with id ", ownerID, " ", ex.what());
          }
      }
  }
  */
  mParts.clear();
  coordinateMap.clear();
}

std::string BondGraph::getId() { return id; };

std::string BondGraph::owner() { return ownerID; }

RCPLIB::RCP<BondGraphInterface> BondGraph::clone() {
  RCPLIB::RCP<BondGraphInterface> bg = createBondGraph();
  RCPLIB::RCP<BondGraph> bgGraph = RCPLIB::rcp_dynamic_cast<BondGraph>(bg);
  bgGraph->cloned = true;
  std::unordered_map<std::string,
                     std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
      componentIdMAP;
  for (auto &c : mComponents) {
    auto copy = createClone(c);
    componentIdMAP[c->getId()] = std::make_tuple(c, copy);
    bgGraph->mComponents.push_back(copy);
  }
  for (auto &b : mBonds) {
    auto fromPort = b->getFromPort();
    auto toPort = b->getToPort();
    auto fromComp = componentIdMAP[fromPort->getComponent()->getId()];
    auto toComp = componentIdMAP[toPort->getComponent()->getId()];
    int opin = 0, ipin = 0;
    // bool oInverting = toPort->getWeight() == -1.0;
    auto oports = std::get<0>(fromComp)->getPorts();
    opin = std::find(oports.begin(), oports.end(), fromPort) - oports.begin();
    auto ports = std::get<0>(toComp)->getPorts();
    ipin = std::find(ports.begin(), ports.end(), toPort) - ports.begin();
    // Creates new bond and returns its handle
    // auto newBond = bgGraph->connect(std::get<1>(fromComp),
    // std::get<1>(toComp), oInverting, opin, ipin);
    auto newBond = bgGraph->connect(std::get<1>(fromComp), std::get<1>(toComp),
                                    opin, ipin);
    newBond->setId(b->getId());
  }
  // ComponentRegistry::defaultRegistry().ownBondgraph(parent,bgGraph);
  return bgGraph;
}

/*! \brief Add a component to the bond graph.
 *  This method add a BondGraphElementBase to the elements list. It will
 *attribute the ID of the passed component based on the current value of the
 *mElementIdCounter and will increment that counter. It will resize the vector
 *if the element ID is greater than the size of the vector. \param inComponent
 *BondGraphElementBase to be added.
 */
void BondGraph::addComponent(const RCPLIB::RCP<BGElement> &inComponent) {
  std::string myID = id;
  if (ownerID != "")
    myID = myID + ":" + ownerID;
  const std::string parentID = inComponent->getParent();
  // Check for circular additions
  if (parentID == myID)
    throw BGException("Component already exists in the bondgraph!");
  std::size_t owned = parentID.find(":");
  if (owned != std::string::npos) {
    std::string ownerid = parentID.substr(owned + 1);
    if (ownerid != "")
      throw BGException(inComponent->getName() + "'s bondgraph (" + parentID +
                        ") is a part of another bondgraph (" + ownerid +
                        "). Parent needs to be free");
  }
  if (parentID == "")
    inComponent->setParent(myID);
  mElementIdCounter++;
  if (inComponent->getId() == "-1") {
    inComponent->setId(mElementIdCounter);
  }
  mComponents.push_back(inComponent);
  ComponentRegistry::defaultRegistry().addOrReplaceComponent(myID, inComponent);
  logDebug("Element ", inComponent->getName(), " of type ",
           inComponent->getType(), " added with parent id ", myID);
}

void BondGraph::addBondGraph(
    const RCPLIB::RCP<BondGraphInterface> &inComponent_) {
  const RCPLIB::RCP<BondGraph> inComponent =
      RCPLIB::rcp_dynamic_cast<BondGraph>(inComponent_);
  // Check if the instance is not owned by anyone else
  if (inComponent->owner() != "")
    throw BGException(
        "Bondgraph is already owned by another instance with id " +
        inComponent->owner());
  /**
   * When adding a bondgraph instance to this instance, proxies for its
   * components are created and added to the current instance. No proxies for
   * bonds are created or new bonds replicating its connection are created. Any
   * changes to the added bond graph are not reflected post this addition
   */

  for (auto &c : inComponent->mComponents) {
    auto copy = createProxy(c); // Create proxies for special handling of
                                // connect
    mComponents.push_back(copy);
  }
  // Proxies and Core components share the same handles to the ports
  // There is no need to create a new bond, use it
  // Bonds may be muted through subsequent connect calls
  for (auto &b : inComponent->mBonds) {
    mBonds.push_back(b);
  }

  inComponent->ownerID = id; // Set the owner
  mParts.push_back(inComponent_);
  *graph += *(inComponent->graph);
  ComponentRegistry::defaultRegistry().addOrReplaceBondgraph(id, inComponent);
}

/**
 * @brief Removes a bondgraph instance with id from itself
 *
 * @param id
 */
void BondGraph::removeBondGraph(std::string id) {
  RCPLIB::RCP<BondGraphInterface> bg;
  bool found = false;
  for (auto c : mParts) {
    if (RCPLIB::rcp_dynamic_cast<BondGraph>(c)->id == id) {
      bg = c;
      found = true;
      break;
    }
  }
  if (found) {
    removeBondGraph(bg);
    auto itr = std::find(mParts.begin(), mParts.end(), bg);
    mParts.erase(itr);
    const RCPLIB::RCP<BondGraph> bgi = RCPLIB::rcp_dynamic_cast<BondGraph>(bg);
    *graph = *graph - *(bgi->graph);
  } else {
    throw BGException("Bondgraph instance with " + id + " not found!");
  }
}

/**
 * @brief Removes a bondgraph instance with RCP handle bg
 *
 * @param bg
 */
void BondGraph::removeBondGraph(const RCPLIB::RCP<BondGraphInterface> &bg_) {
  const RCPLIB::RCP<BondGraph> bg = RCPLIB::rcp_dynamic_cast<BondGraph>(bg_);
  std::map<std::string, std::vector<RCPLIB::RCP<BGElement>>::iterator>
      rComponents; //!< List of all elements
  for (std::vector<RCPLIB::RCP<BGElement>>::iterator c = mComponents.begin();
       c != mComponents.end(); c++) {
    rComponents[(*c)->getId()] = c;
  }
  std::vector<std::vector<RCPLIB::RCP<BGElement>>::iterator> items;
  for (auto c : bg->mComponents) {
    auto res = rComponents.find(c->getId());
    if (res != rComponents.end()) {
      items.push_back(rComponents[c->getId()]);
    }
  }
  std::vector<std::vector<RCPLIB::RCP<BondInterface>>::iterator> bonds;
  for (std::vector<RCPLIB::RCP<BondInterface>>::iterator c = bg->mBonds.begin();
       c != bg->mBonds.end(); c++) {
    auto res = std::find(mBonds.begin(), mBonds.end(), *c);
    if (res != mBonds.end()) {
      bonds.push_back(res);
    }
  }
  for (auto &itr : items) {
    mComponents.erase(itr);
  }
  for (auto &bir : bonds) {
    mBonds.erase(bir);
  }
}

/*! \brief Add a bond to the bond graph.
 *  This method add a Bond to the bond list. It will attribute the ID
 *	of the passed bond based on the current value of the mBondIdCounter
 *	and will then increment that counter. It will resize the vector if the
 *bond ID is greater than the size of the vector. \param  inBond Bond to be
 *added
 */
void BondGraph::addBond(const RCPLIB::RCP<BondInterface> &inBond) {
  if (inBond->getId() < 0) {
    inBond->setId(mBondIdCounter++);
  }

  mBonds.push_back(inBond);
}

/*! \brief Connect two components.
 *  This method define the connectivity between two BondGraphElementBase. The
 *power direction is from inOrigin to inDest. The causality is the one defined
 *at the destination component. Therefore, if the causality is said to be
 *eEffortCausal, the causality stroke will be applied to the inDest side of the
 *bond. Using this method, a new bond is created. \param  inOrigin
 *BondGraphElementBase at the origin of the bond \param  inDest
 *BondGraphElementBase at the end of the bond \param  inCausality Causality
 *relation between the component, looking at inDest causality. \param	opin
 *Select the port for origin connection. \param	ipin       Select the port for
 *destination connection. \return	The newly created bond.
 */
RCPLIB::RCP<BondInterface>
BondGraph::connect(const RCPLIB::RCP<BGElement> &inOrigin_,
                   const RCPLIB::RCP<BGElement> &inDest_, int opin_,
                   int ipin_) {
  RCPLIB::RCP<BondGraphElementBase> inOrigin =
      RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(inOrigin_);
  RCPLIB::RCP<BondGraphElementBase> inDest =
      RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(inDest_);
  int opin = opin_;
  int ipin = ipin_;
  // Create both ports
  RCPLIB::RCP<PortInterface> lFromPort, lToPort;
  auto &oPorts = inOrigin->getPorts();
  auto &iPorts = inDest->getPorts();

  // Replace non junction elements, single port elements that do not belong to
  // this Bondgraph and connect to their Junction For reactions connect to the
  // port
  if (inOrigin->getParent() != "" && inOrigin->getComponentGroup() != eJ &&
      inOrigin->getParent() != id) {
    std::ostringstream ss;
    if (oPorts.size() < opin) {
      ss << "Incorrect portnumber for origin element " << inOrigin->getName()
         << " Parent bondgraph " << inOrigin->getParent();
      std::string msg = ss.str();
      throw BGException(msg);
    }
    const RCPLIB::RCP<BondInterface> &bond = oPorts[opin]->getBond();
    if (oPorts.size() ==
        1) { // Single port element, find the connecting element from the bond
      const RCPLIB::RCP<PortInterface> &rp = oPorts[opin];
      const RCPLIB::RCP<PortInterface> &oend = bond->getOtherEnd(rp);
      inOrigin =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(oend->getComponent());
      ss << inOrigin_->getName() << " component released, " << inDest->getName()
         << " will be connected to " << inOrigin->getName();
      oPorts = inOrigin->getPorts();
      // Ensure that port number matches to that of the one being released
      if (oPorts.size() > 1) {
        opin = std::find(oPorts.begin(), oPorts.end(), rp) - oPorts.begin();
        ss << " at port number " << opin;
      }
      std::string msg = ss.str();
      logInfo(msg);
    }
    removeBond(bond); // This will disconnect the ports and components
  }
  if (inDest->getParent() != "" && inDest->getComponentGroup() != eJ &&
      inDest->getParent() != id) {
    std::ostringstream ss;
    if (iPorts.size() < ipin) {
      ss << "Incorrect portnumber for destination element " << inDest->getName()
         << " Parent bondgraph " << inDest->getParent();
      std::string msg = ss.str();
      throw BGException(msg);
    }
    const RCPLIB::RCP<BondInterface> &bond = iPorts[ipin]->getBond();
    if (oPorts.size() ==
        1) { // Single port element, find the connecting element from the bond
      const RCPLIB::RCP<PortInterface> rp = iPorts[ipin];
      const RCPLIB::RCP<PortInterface> &oend = bond->getOtherEnd(rp);
      inDest =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(oend->getComponent());
      ss << inDest_->getName() + " component released, " + inOrigin->getName() +
                " will be connected to " + inDest->getName();
      iPorts = inDest->getPorts();
      // Ensure that port number matches to that of the one being released
      if (oPorts.size() > 1) {
        ipin = std::find(iPorts.begin(), iPorts.end(), rp) - iPorts.begin();
        ss << " at port number " << ipin;
      }
      std::string msg = ss.str();
      logInfo(msg);
    }
    removeBond(bond); // This will disconnect the ports and components
  }

  // Create a new bond
  const RCPLIB::RCP<BondInterface> &lBond = createBond(-1);
  addBond(lBond);

  if (inOrigin->getComponentGroup() != eJ) {
    if (oPorts.size() < opin) {
      std::ostringstream ss;
      ss << "Incorrect portnumber for origin element " << inOrigin->getName();
      std::string msg = ss.str();
      throw BGException(msg);
    }
    if (oPorts.size() == 0 && opin == 0) {
      lFromPort = createPort(false);
      inOrigin->connect(lFromPort); // Add to element port list, which then
                                    // links the element with the port
    } else {
      lFromPort = oPorts[opin];
      lFromPort->connect(
          inOrigin
              .create_weak()); // Link element with port, passing a strong
                               // pointer leads to inOrigin to not be released
    }
  } else {
    lFromPort = createPort(false);
    inOrigin->connect(lFromPort, opin);
  }
  if (inDest->getComponentGroup() != eJ) {
    if (iPorts.size() < ipin) {
      std::ostringstream ss;
      ss << "Incorrect portnumber for destination element "
         << inDest->getName();
      std::string msg = ss.str();
      throw BGException(msg);
    }
    if (iPorts.size() == 0 && ipin == 0) {
      lToPort = createPort(true);
      inDest->connect(lToPort); // Add to element port list, which then links
                                // the element with the port
    } else {
      lToPort = iPorts[ipin];
      lToPort->connect(
          inDest.create_weak()); // Link element with port, passing a strong
                                 // pointer leads to inOrigin to not be released
    }
  } else {
    lToPort = createPort(true);
    inDest->connect(lToPort, ipin); // Add to element port list, which then
                                    // links the element with the port
  }

  if (inOrigin->getType() == BG::PassiveType::eOne) {
    lFromPort->setWeight(-1.0);
  }

  // Connect the elements together
  lBond->connect(lFromPort, lToPort);
  graph->insert_edge(lFromPort->getComponent()->getId(),
                     lToPort->getComponent()->getId());
  return lBond;
}

void BondGraph::removeBond(const RCPLIB::RCP<BondInterface> &inBond) {
  inBond->getFromPort()->getComponent()->disconnect(inBond->getFromPort());
  inBond->getToPort()->getComponent()->disconnect(inBond->getToPort());

  // Delete bond
  // Didnt work as Teuchos pointer instances are different although internal
  // pointers are same
  /*
  auto lBondIter = std::find(mBonds.begin(), mBonds.end(), inBond);

  if (lBondIter != mBonds.end()) {
      graph.remove_edge((*lBondIter)->getFromPort()->getComponent()->getId(),
  (*lBondIter)->getToPort()->getComponent()->getId()); mBonds.erase(lBondIter);
  }
  */

  // std::vector<RCPLIB::RCP<BondInterface>>::iterator lBondIter =
  // mBonds.begin();
  for (auto lBondIter = mBonds.begin(); lBondIter != mBonds.end();
       lBondIter++) {
    const RCPLIB::RCP<Bond> bi = RCPLIB::rcp_dynamic_cast<Bond>(*lBondIter);
    const RCPLIB::RCP<Bond> ib = RCPLIB::rcp_dynamic_cast<Bond>(inBond);
    if (*bi == *ib) {
      graph->remove_edge((inBond)->getFromPort()->getComponent()->getId(),
                         (inBond)->getToPort()->getComponent()->getId());
      mBonds.erase(lBondIter);
      break;
    }
  }
  // if (lBondIter != mBonds.end()) {
  //   mBonds.erase(lBondIter);
  // }
}

void BondGraph::removeComponent(
    std::vector<RCPLIB::RCP<BGElement>> &inComponents) {
  for (unsigned int i = 0; i < inComponents.size(); ++i) {
    removeComponent(inComponents[i]);
  }
}

void BondGraph::removeComponent(const RCPLIB::RCP<BGElement> &inComponent) {
  std::vector<RCPLIB::RCP<BGElement>>::iterator lComponentIter =
      find(mComponents.begin(), mComponents.end(), inComponent);
  assert(lComponentIter != mComponents.end());

  for (unsigned int i = 0; i < inComponent->getPorts().size(); ++i) {
    const RCPLIB::RCP<PortInterface> &lPort = inComponent->getPorts()[i];
    const RCPLIB::RCP<BondInterface> &lBond = lPort->getBond();

    // Delete bond
    removeBond(lBond);
  }
  mComponents.erase(lComponentIter);
}

void BondGraph::reconnectComponent(
    const RCPLIB::RCP<PortInterface> &inComponentPort,
    const RCPLIB::RCP<BGElement> &inNewOrigin) {
  const RCPLIB::RCP<BondInterface> &lBond = inComponentPort->getBond();
  graph->remove_edge(lBond->getFromPort()->getComponent()->getId(),
                     lBond->getToPort()->getComponent()->getId());

  // Remove connection to the old origin component
  const RCPLIB::RCP<PortInterface> &lMovingPort =
      lBond->getOtherEnd(inComponentPort);
  lMovingPort->getComponent()->disconnect(lMovingPort);
  // Add a new connection to the new origin component
  inNewOrigin->connect(lMovingPort);
  graph->insert_edge(lBond->getFromPort()->getComponent()->getId(),
                     lBond->getToPort()->getComponent()->getId());
}

std::string BondGraph::serialize(bool inIndent) const { return ""; }

std::string BondGraph::serializeAsCellML(std::string modelName) {
  auto result = computeStateEquation();
  // Filename, contents
  const RCPLIB::RCP<BondGraphInterface> interface = RCPLIB::rcp(this);
  std::map<std::string, std::string> res =
      getCellML(modelName, interface, result);
  return "";
}

void BondGraph::read(std::string filename) {}

RCPLIB::RCP<BondInterface>
BondGraph::connect(const RCPLIB::RCP<BGElement> &inOrigin, int opin,
                   const RCPLIB::RCP<BGElement> &inDest, int ipin) {
  return connect(inOrigin, inDest, opin, ipin);
};
RCPLIB::RCP<BondInterface>
BondGraph::connect(const RCPLIB::RCP<BGElement> &inOrigin, int opin,
                   const RCPLIB::RCP<BGElement> &inDest) {
  return connect(inOrigin, inDest, opin, 0);
};
RCPLIB::RCP<BondInterface>
BondGraph::connect(const RCPLIB::RCP<BGElement> &inOrigin,
                   const RCPLIB::RCP<BGElement> &inDest, int ipin) {
  return connect(inOrigin, inDest, 0, ipin);
};
RCPLIB::RCP<BondInterface>
BondGraph::connectInverting(const RCPLIB::RCP<BGElement> &inOrigin, int opin,
                            const RCPLIB::RCP<BGElement> &inDest) {
  return connect(inOrigin, inDest, opin, 0);
};
RCPLIB::RCP<BondInterface>
BondGraph::connectInverting(const RCPLIB::RCP<BGElement> &inOrigin,
                            const RCPLIB::RCP<BGElement> &inDest, int ipin) {
  return connect(inOrigin, inDest, 0, ipin);
};

std::vector<std::tuple<std::string, const RCPLIB::RCP<BGElement> &, bool>>
BondGraph::getComponents() {
  std::vector<std::tuple<std::string, const RCPLIB::RCP<BGElement> &, bool>>
      result;
  for (auto c : mComponents) {
    std::string readableName =
        ComponentRegistry::defaultRegistry().nameMap[c->getId()];
    if (readableName.size() == 0) {
      // Set it to the element name - when user sets name
      readableName = c->getId();
    }
    bool proxy = c->isProxy();
    // Check for passing c as a local variable and loss of memory mapping
    result.push_back(std::tie(readableName, c, proxy));
  }
  return result;
}

const RCPLIB::RCP<BGElement> &BondGraph::getElementByName(std::string name) {
  for (int i = 0; i < mComponents.size(); i++) {
    if (mComponents[i]->getName() == name) {
      return mComponents[i];
    }
  }
  throw BGException("Element with name " + name + " not found");
}

std::vector<RCPLIB::RCP<BGElement>> &BondGraph::getSources() {
  return mSources;
}
const std::vector<RCPLIB::RCP<BondInterface>> &BondGraph::getBonds() const {
  return mBonds;
}

RCPLIB::RCP<BondGraphInterface> createBondGraph() {
  RCPLIB::RCP<BondGraphInterface> ptr = RCPLIB::rcp(new BondGraph());
  ComponentRegistry::defaultRegistry().ownBondgraph("", ptr);
  return ptr;
}

RCPLIB::RCP<BGElement> __createResistance_do_not_call_outside_of_lib(
    const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "R_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Resistance(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement> __createCapacitance_do_not_call_outside_of_lib(
    const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "C_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Capacitance(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

#include "factorymethods.cpp"

RCPLIB::RCP<BGElement> __createInductance_do_not_call_outside_of_lib(
    const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "I_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Inductance(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createPotentialSource(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Se_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new PotentialSource(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createFlowSource(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Sf_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new FlowSource(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createConcentration(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Ce_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Concentration(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createTransformer(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Tr_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Transformer(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement> createGyrator(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Gy_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Gyrator(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement> createReaction(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Rx_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new Reaction(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createChemostat(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Ss_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new ChemicalSource(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement> createFlowstat(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Sf_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new ChemicalFlux(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createUserDefined(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "Ue_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new UserDefined(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createPortHamiltonian(std::string expr, std::vector<std::string> states,
                      const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "PH_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  // std::vector<std::string> states_(states);
  if (proxy.is_null()) {
    if (states.size() == 0 || expr == "") {
      throw BGException(
          "PortHamiltonion's need an energy function and state description");
    }
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
    // This is not required as PortHamiltonian will return
    // after the copying of elementImpl
    /*
    for(int i=0;i<data->numStates;i++){
        std::string& sp = mParameter[i]->prefix;
        std::string sn = sp.substr(0,sp.rfind("_"));
        states_.push_back(sn);
    }*/
  }
  // If creating a proxy, the constructor will return post copying the
  // elementImpl - no processing is done so states and expr can be empty
  auto ptr =
      RCPLIB::rcp(new PortHamiltonian(expr, states, data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createOneJunction(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "1J_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new OneJunction(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement>
createZeroJunction(const RCPLIB::RCP<BGElementData> &proxy) {
  static unsigned int ctr = 0;
  std::string readableName = "0J_" + std::to_string(++ctr);
  RCPLIB::RCP<BGElementData> data = proxy;
  std::string id = "-1";
  bool proxyFlag = true;
  if (proxy.is_null()) {
    data = RCPLIB::rcp(new BGElementData);
    proxyFlag = false;
  } else {
    id = data->mId;
  }
  auto ptr = RCPLIB::rcp(new ZeroJunction(data, id, proxyFlag));
  ptr->rcpPtr = ptr.create_weak();
  if (!proxyFlag)
    ComponentRegistry::defaultRegistry().ownComponent("", ptr, readableName);
  return ptr;
}

RCPLIB::RCP<BGElement> createProxy(const RCPLIB::RCP<BGElement> &obj_) {
  // auto ptr = RCPLIB::rcp(new BondGraphElementBase(obj));
  const RCPLIB::RCP<BondGraphElementBase> obj =
      RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(obj_);
  RCPLIB::RCP<BondGraphElementBase> ptr;
  auto id = obj->mId;
  std::string parent = obj->parent;
  const RCPLIB::RCP<BGElementData> data = obj->_data;
  PassiveType etype = obj->getType();
// Define
#include "createElements.func"

  ptr->mId = id;
  ptr->parent = parent;
  logDebug("Created a proxy for element ", data->mName, " type ",
           data->mElementType, " parent ", parent);
  return ptr;
}

RCPLIB::RCP<BGElement> createClone(const RCPLIB::RCP<BGElement> &obj_) {
  const RCPLIB::RCP<BondGraphElementBase> obj =
      RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(obj_);
  RCPLIB::RCP<BondGraphElementBase> ptr;
  const RCPLIB::RCP<BGElementData> data = RCPLIB::null;
  PassiveType etype = obj->getType();

// Define
#include "createElements.func"

  // Update parameters, for port hamiltonians, not required as they are computed
  if (obj->getType() != ePortHamiltonian) {
    const RCPLIB::RCP<BGElementData> data = ptr->_data;
    data->mName = obj->mName;
    data->mElementType = obj->mElementType;
    data->mComponentGroup = obj->mComponentGroup;
    data->mDomain = obj->mDomain;
    data->numStates = obj->numStates;
    data->hamiltonian = obj->hamiltonian;
    data->allocatedPorts = obj->allocatedPorts;
    data->constitutiveEq = std::vector<std::string>(obj->constitutiveEq);
    data->constitutiveEqIndex =
        std::vector<int>(obj->constitutiveEq.size(), -1);

    for (int c = 0; c < data->mParameter.size(); c++) {
      data->mParameter[c] = obj->mParameter[c];
    }
    logDebug("Created a clone for element ", data->mName, " type ",
             data->mElementType);
  } else {
    // Note data in global context is null
    logDebug("Created a clone for element ", ptr->mName, " type ",
             ptr->mElementType);
  }
  return ptr;
}

RCPLIB::RCP<BGElement> createBondgraphElement(std::string methodName,
                                              std::string expr) {
  CALL_FACTORY_METHODS_BY_NAME
  throw BGException("FactoryMethod " + methodName + " not found");
}

#include "computeStateEquation.h"

#include "computePortHamiltonian.h"

std::string checkUnits(std::string unit, std::string target) {
  units::precise_unit compare = units::unit_from_string(unit);
  units::precise_unit targetUnit = units::unit_from_string(target);
  if (compare.base_units() !=
      targetUnit.base_units()) { // Compare the base units as multiplies can be
                                 // different
    return units::to_string(compare);
  }
  return "";
}

} // namespace BG
