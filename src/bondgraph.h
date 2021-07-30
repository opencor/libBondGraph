#ifndef BONDGRAPH__H_
#define BONDGRAPH__H_
#include "Bond.h"
#include "Port.h"
#include "RCP.h"
#include "export.h"
#include <symengine/matrix.h>
#include <tuple>
#include <vector>
#include "thirdparty/ngraph.hpp"

namespace BG {

/*! \brief %Bond Graph
		This class is a collection of elements and bonds. 
		To define a bond graph, individual elements need to be instantiated and added to the collection 
		The elements should then be connected together using the connect method.

		Here is an example how to create a bond graph:
		\code
		//Create the source
        auto ioBondGraph = createBondGraph();

        //Create the storage
        auto A = createConcentration();
        ioBondGraph->addComponent(A);

        auto B = createConcentration();
        ioBondGraph->addComponent(B);
        
        //Create the reaction
        auto Re = createReaction();
        ioBondGraph->addComponent(Re);

        //Create the junctions
        auto Y_A = createOneJunction();
        ioBondGraph->addComponent(Y_A);

        auto Y_B = createOneJunction();
        ioBondGraph->addComponent(Y_B);

        //Create the bonds
        ioBondGraph->connect(A, Y_A);
        ioBondGraph->connect(B, Y_B);
        ioBondGraph->connectInverting(Re, 0, Y_A);
        ioBondGraph->connectInverting(Re, 1, Y_B);
		\endcode

		The bond graph is represented as a linked list. All elements of the bond
		graph are interconnected through ports. Each port connects to a Component via a Bond.
		Therefore, two component are linked like this : 
			%Element -> %Port -> %Bond -> %Port -> %Element
	*/

class EXPORTED BondGraph
{
private:
    std::string id; //! Unique identifier for the bondgraph
    bool cloned = false; //! True when bondgraph is cloned, getComponents will return proxies
    std::string ownerID = ""; //! ID of the bondgraph instance, if this instance is a component in another bondgraph with id=ownerID, sanity checking to avoid sharing one instance with multiple parents
protected:
    unsigned int mElementIdCounter = 1; //!< Id counter for elements
    unsigned int mBondIdCounter = 1; //!< Id counter for bonds
    NGraph::sGraph graph; //Graph data structure to capture connectivity
    std::vector<RCPLIB::RCP<BondGraphElementBase>> mComponents; //!< List of all elements
    std::vector<RCPLIB::RCP<BondGraphElementBase>> mSources; //!< List of all sources
    std::vector<RCPLIB::RCP<Bond> > mBonds; //!< List of all bonds
    std::vector<RCPLIB::RCP<BondGraph>> mParts;

    //State equation related members and methods
    int mUcount = -1; //!< Number of source bond
    int mScount = -1; //!< Number of storage bond with integral causality
    int mSDcount = -1; //!< Number of storage bond with derivative causality
    int mRcount = -1; //!< Number of dissipative bond
    int mJcount = -1; //!< Number of junction bond

    //Simulation related members
    std::vector<double> getSourceValues(double inTime) const;
    //Causality assignment related members and methods
    void countComponentType();
    // Coordinates
    SymEngine::DenseMatrix computedCoordinates;
    std::unordered_map<std::string, long int> coordinateMap;

    //Bondgraphs can only be created by the factory method
    BondGraph();
    //Support function for bond creation
    void addBond(const RCPLIB::RCP<Bond> &inBond);
    BondGraph(const BondGraph &inBondGraph) = delete;
    void operator=(const BondGraph &) = delete;

public:
    virtual ~BondGraph();

    std::string getId(){
        return id;
    };

    /**
     * @brief ID of the bondgraph that owns this bondgraph instance, empty if its independent
     * 
     * @return std::string 
     */
    std::string owner();
    /**
     * @brief Create a clone of the input bondgraph
     * A copy of the bondgraph is made by copying all elements (proxies are also concretely allocated) and bonds are restablised 
     * 
     * @return RCPLIB::RCP<BondGraph> 
     */
    RCPLIB::RCP<BondGraph> clone();
    //Creation related methods
    /**
     * @brief Add a bondgraph element to the graph
     * 
     * @param inComponent Element to be added
     */
    void addComponent(const RCPLIB::RCP<BondGraphElementBase> &inComponent);
    //Creation related methods
    /**
     * @brief Add a bondgraph to the graph
     * 
     * @param inComponent Element to be added
     */
    void addBondGraph(const RCPLIB::RCP<BondGraph> &inComponent);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param inDest   The element to which power is expected to flow in
     * @param oInverting When one of the elements has two ports (input and output), the other element is connected to the input port of this element of oInverting is true. default is false 
     * @param opin Specify the port from which power should flow out: default 0
     * @param ipin Specify the port to  which power should flow in: default 0
     * @return RCPLIB::RCP<Bond> The bond connecting the two elements
     */
    RCPLIB::RCP<Bond> connect(const RCPLIB::RCP<BondGraphElementBase> &inOrigin, const RCPLIB::RCP<BondGraphElementBase> &inDest, bool oInverting = false, int opin = 0, int ipin = 0);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param opin Specify the port from which power should flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin Specify the port to  which power should flow in
     * @param oInverting When one of the elements has two ports (input and output), the other element is connected to the input port of this element of oInverting is true
     * @return RCPLIB::RCP<Bond> The bond connecting the two elements
     */
    RCPLIB::RCP<Bond> connect(const RCPLIB::RCP<BondGraphElementBase> &inOrigin, int opin, const RCPLIB::RCP<BondGraphElementBase> &inDest, int ipin, bool oInverting);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param opin Specify the port from which power should flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin Specify the port to  which power should flow in: default 0
     * @param oInverting When one of the elements has two ports (input and output), the other element is connected to the input port of this element of oInverting is true. default is false 
     * @return RCPLIB::RCP<Bond> The bond connecting the two elements
     */
    RCPLIB::RCP<Bond> connect(const RCPLIB::RCP<BondGraphElementBase> &inOrigin, int opin, const RCPLIB::RCP<BondGraphElementBase> &inDest, bool oInverting = false);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin Specify the port to  which power should flow in
     * @param oInverting When one of the elements has two ports (input and output), the other element is connected to the input port of this element of oInverting is true. default is false      
     * @return RCPLIB::RCP<Bond> The bond connecting the two elements
     */
    RCPLIB::RCP<Bond> connect(const RCPLIB::RCP<BondGraphElementBase> &inOrigin, const RCPLIB::RCP<BondGraphElementBase> &inDest, int ipin, bool oInverting = false);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param opin     Specify the port from which power should flow out
     * @param inDest   The element to which power is expected to flow in
     * @return RCPLIB::RCP<Bond> The bond connecting the two elements
     */
    RCPLIB::RCP<Bond> connectInverting(const RCPLIB::RCP<BondGraphElementBase> &inOrigin, int opin, const RCPLIB::RCP<BondGraphElementBase> &inDest);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin     Specify the port from which power should flow in
     * @return RCPLIB::RCP<Bond> The bond connecting the two elements
     */
    RCPLIB::RCP<Bond> connectInverting(const RCPLIB::RCP<BondGraphElementBase> &inOrigin, const RCPLIB::RCP<BondGraphElementBase> &inDest, int ipin);

    /**
     * @brief Removes a bondgraph instance with id from itself
     * 
     * @param id 
     */
    void removeBondGraph(std::string id);

    /**
     * @brief Removes a bondgraph instance with RCP handle bg
     * 
     * @param bg 
     */
    void removeBondGraph(const RCPLIB::RCP<BondGraph> &bg);

    //Bondgraph manipulation methods
    /**
     * @brief Remove a bond
     * 
     * @param inBond handle to bond
     */
    virtual void removeBond(const RCPLIB::RCP<Bond> &inBond);
    /**
     * @brief Remove component and associated linkages
     * 
     * @param inComponent component to be removed
     */
    virtual void removeComponent(const RCPLIB::RCP<BondGraphElementBase> &inComponent);
    /**
     * @brief Remove all components and their associated linkages in the list
     * 
     * @param inComponents std::vector of element handles that need to be removed
     */
    virtual void removeComponent(std::vector<RCPLIB::RCP<BondGraphElementBase>> &inComponents);
    /**
     * @brief Reconnect a component to a port
     * 
     * @param inComponentPort Port to which the component needs be connected
     * @param inNewOrigin     Component that needs to be connected
     */
    virtual void reconnectComponent(const RCPLIB::RCP<Port> &inComponentPort, const RCPLIB::RCP<BondGraphElementBase> &inNewOrigin);

    //State equation related methods
    /**
     * @brief Compute the state equations for the bondgraph 
     * First argument of the tuple will be false if the logic concludes that the bg is incorrect, for instance, some state and power variables are evaluating to zero, it will still return the equations
     * Second argument: system of equations corresponding to each dof - dof_symbol, eqn
     * Third argument: system of equations corresponding to each bond - bond_variable, eqn
     * Fourth argument: list of constraints
     * Firth argument: physical dimensions, value and variable type associated with every the model: the last element identifies if it is a state (s), control (c) or parameter (p)
     * @return std::tuple<bool,SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::string>> 
     */
    std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>> > computeStateEquation();
    /**
     * @brief Get the Components object
     * Returns a tuple with either the readableName (if available) or id, handle to it and whether it is a proxy
     * @return std::vector<std::tuple<str::string,RCPLIB::RCP<BondGraphElementBase>,bool>>& 
     */
    std::vector<std::tuple<std::string,const RCPLIB::RCP<BondGraphElementBase>&,bool>> getComponents();
    /**
     * @brief Get the list of source elements associated with this bondgraph, entries are valid only after computeStateEquation is called
     * 
     * @return std::vector<RCPLIB::RCP<BondGraphElementBase>>& List of elements
     */
    virtual std::vector<RCPLIB::RCP<BondGraphElementBase>> &getSources();
    /**
     * @brief Get the list of Bonds associated with this bondgraph, entries are valid only after computeStateEquation is called
     * 
     * @return const std::vector<RCPLIB::RCP<Bond>>& List of bonds
     */
    virtual const std::vector<RCPLIB::RCP<Bond>> &getBonds() const;

    //IO related methods
    //!Serialize the bond graph in JSON form.
    std::string serialize(bool inIndent = false) const;
    //!Serialize the bond graph in flattend one component CellML form.
    std::string serializeAsCellML(std::string modelName);
    virtual void read(std::string filename);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraph> &p);
    friend void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraph> &p);
    friend RCPLIB::RCP<BondGraph> createBondGraph();
};

/**
 * @brief Create a Bond Graph instance
 * 
 * @return Reference counted pointer to the bondgraph 
 */
EXPORTED RCPLIB::RCP<BondGraph> createBondGraph();

/**
 * @brief Create a Resistane instance
 * Library users are not expected to call this as element definition is incomplete and domain is not assigned.
 * Users should call one of the generated factory methods
 * @return Reference counted pointer to the element 
 */
RCPLIB::RCP<Resistance> __createResistance_do_not_call_outside_of_lib(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Capacitance instance
 * Library users are not expected to call this as element definition is incomplete and domain is not assigned.
 * Users should call one of the generated factory methods 
 * @return Reference counted pointer to the element 
 */
RCPLIB::RCP<Capacitance> __createCapacitance_do_not_call_outside_of_lib(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Inductance instance
 * Library users are not expected to call this as element definition is incomplete and domain is not assigned.
 * Users should call one of the generated factory methods
 * @return Reference counted pointer to the element 
 */
RCPLIB::RCP<Inductance> __createInductance_do_not_call_outside_of_lib(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

#include "factorymethods.h"

/**
 * @brief Create a Effort Source instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<PotentialSource> createPotentialSource(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Flow Source instance
 * 
 * @return Reference counted pointer to the element 
 */
EXPORTED RCPLIB::RCP<FlowSource> createFlowSource(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Concentration instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<Concentration> createConcentration(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Transformer instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<Transformer> createTransformer(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Gyrator instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<Gyrator> createGyrator(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Reaction instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<Reaction> createReaction(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a User defined instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<UserDefined> createUserDefined(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a User defined instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<PortHamiltonian> createPortHamiltonian(std::string expr="",
                                                            std::vector<std::string> states={}, const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a One Junction instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<OneJunction> createOneJunction(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create a Zero Junction instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<ZeroJunction> createZeroJunction(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

/**
 * @brief Create proxy instance
 * 
 * @return Reference counted pointer to the element
 */
EXPORTED RCPLIB::RCP<BondGraphElementBase> createProxy(const RCPLIB::RCP<BondGraphElementBase> &obj);

/**
 * @brief Create a Clone object
 * 
 * @param Reference counted pointer to the element
 * @return RCPLIB::RCP<BondGraphElementBase> 
 */
EXPORTED RCPLIB::RCP<BondGraphElementBase> createClone(const RCPLIB::RCP<BondGraphElementBase> &obj);

/**
 * @brief Create a Bondgraph element by calling the factory method specified by name
 * 
 * @param Reference counted pointer to the element
 * @return RCPLIB::RCP<BondGraphElementBase> 
 */
EXPORTED RCPLIB::RCP<BondGraphElementBase> createBondgraphElement(std::string methodName);
/**
 * @brief Clears all bondgraph(s), annotation(s) and elements from the registry and frees associated memory
 *  User is expected to save any changes to these instances prior to this call
 */
EXPORTED void newWorkSpace();

} // namespace BG

#endif