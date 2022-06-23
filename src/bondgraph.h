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

#include "bondgraph.hpp"
#include "Bond.h"
#include "Port.h"
#include "RCP.h"
#include "export.h"
#include <tuple>
#include <vector>


namespace NGraph{
    template <typename T>
    class tGraph;
    typedef tGraph<std::string> sGraph;
}




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

class BondGraph : public BondGraphInterface
{
private:
    std::string id; //! Unique identifier for the bondgraph
    bool cloned = false; //! True when bondgraph is cloned, getComponents will return proxies
    std::string ownerID = ""; //! ID of the bondgraph instance, if this instance is a component in another bondgraph with id=ownerID, sanity checking to avoid sharing one instance with multiple parents
protected:
    unsigned int mElementIdCounter = 1; //!< Id counter for elements
    unsigned int mBondIdCounter = 1; //!< Id counter for bonds
    RCPLIB::RCP<NGraph::sGraph> graph; //Graph data structure to capture connectivity
    std::vector<RCPLIB::RCP<BGElement> > mComponents; //!< List of all elements
    std::vector<RCPLIB::RCP<BGElement> > mSources; //!< List of all sources
    std::vector<RCPLIB::RCP<BondInterface> > mBonds; //!< List of all bonds
    std::vector<RCPLIB::RCP<BondGraphInterface>> mParts;

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
    RCPLIB::RCP<SymEngine::DenseMatrix> computedCoordinates;
    std::unordered_map<std::string, long int> coordinateMap;

    //Bondgraphs can only be created by the factory method
    BondGraph();
    //Support function for bond creation
    void addBond(const RCPLIB::RCP<BondInterface> &inBond);
    BondGraph(const BondGraph &inBondGraph) = delete;
    void operator=(const BondGraph &) = delete;

public:
    virtual ~BondGraph();
    /**
     * @brief Get the unique Id object
     * 
     * @return std::string 
     */
     virtual std::string getId();
    /**
     * @brief ID of the bondgraph that owns this bondgraph instance, empty if its independent
     * 
     * @return std::string 
     */
     virtual std::string owner();
    /**
     * @brief Create a clone of the input bondgraph
     * A copy of the bondgraph is made by copying all elements (proxies are also concretely allocated) and bonds are restablised 
     * 
     * @return RCPLIB::RCP<BondGraphInterface> 
     */
     virtual RCPLIB::RCP<BondGraphInterface> clone();
    //Creation related methods
    /**
     * @brief Add a bondgraph element to the graph
     * 
     * @param inComponent Element to be added
     */
     virtual void addComponent(const RCPLIB::RCP<BGElement>  &inComponent);
    //Creation related methods
    /**
     * @brief Add a bondgraph to the graph
     * 
     * @param inComponent Element to be added
     */
     virtual void addBondGraph(const RCPLIB::RCP<BondGraphInterface> &inComponent);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param inDest   The element to which power is expected to flow in
     * @param opin Specify the port from which power should flow out: default 0
     * @param ipin Specify the port to  which power should flow in: default 0
     * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
     */
     virtual RCPLIB::RCP<BondInterface> connect(const RCPLIB::RCP<BGElement>  &inOrigin, const RCPLIB::RCP<BGElement>  &inDest, int opin = 0, int ipin = 0);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param opin Specify the port from which power should flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin Specify the port to  which power should flow in
     * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
     */
     virtual RCPLIB::RCP<BondInterface> connect(const RCPLIB::RCP<BGElement>  &inOrigin, int opin, const RCPLIB::RCP<BGElement>  &inDest, int ipin);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param opin Specify the port from which power should flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin Specify the port to  which power should flow in: default 0
     * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
     */
     virtual RCPLIB::RCP<BondInterface> connect(const RCPLIB::RCP<BGElement>  &inOrigin, int opin, const RCPLIB::RCP<BGElement>  &inDest);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin Specify the port to  which power should flow in
     * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
     */
     virtual RCPLIB::RCP<BondInterface> connect(const RCPLIB::RCP<BGElement>  &inOrigin, const RCPLIB::RCP<BGElement>  &inDest, int ipin);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param opin     Specify the port from which power should flow out
     * @param inDest   The element to which power is expected to flow in
     * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
     */
     virtual RCPLIB::RCP<BondInterface> connectInverting(const RCPLIB::RCP<BGElement>  &inOrigin, int opin, const RCPLIB::RCP<BGElement>  &inDest);
    /**
     * @brief Connect the power ports of two bondgraph elements that exist within the bondgraph
     * 
     * @param inOrigin The element from which power is expected to flow out
     * @param inDest   The element to which power is expected to flow in
     * @param ipin     Specify the port from which power should flow in
     * @return RCPLIB::RCP<BondInterface> The bond connecting the two elements
     */
     virtual RCPLIB::RCP<BondInterface> connectInverting(const RCPLIB::RCP<BGElement>  &inOrigin, const RCPLIB::RCP<BGElement>  &inDest, int ipin);

    /**
     * @brief Removes a bondgraph instance with id from itself
     * 
     * @param id 
     */
     virtual void removeBondGraph(std::string id);

    /**
     * @brief Removes a bondgraph instance with RCP handle bg
     * 
     * @param bg 
     */
     virtual void removeBondGraph(const RCPLIB::RCP<BondGraphInterface> &bg);

    //Bondgraph manipulation methods
    /**
     * @brief Remove a bond
     * 
     * @param inBond handle to bond
     */
     virtual void removeBond(const RCPLIB::RCP<BondInterface> &inBond);
    /**
     * @brief Remove component and associated linkages
     * 
     * @param inComponent component to be removed
     */
     virtual void removeComponent(const RCPLIB::RCP<BGElement>  &inComponent);
    /**
     * @brief Remove all components and their associated linkages in the list
     * 
     * @param inComponents std::vector of element handles that need to be removed
     */
     virtual void removeComponent(std::vector<RCPLIB::RCP<BGElement> > &inComponents);
    /**
     * @brief Reconnect a component to a port
     * 
     * @param inComponentPort Port to which the component needs be connected
     * @param inNewOrigin     Component that needs to be connected
     */
     virtual void reconnectComponent(const RCPLIB::RCP<PortInterface> &inComponentPort, const RCPLIB::RCP<BGElement>  &inNewOrigin);

    //State equation related methods
    /**
     * @brief Compute the state equations for the bondgraph 
     * @return instance of ComputeEquationResults
     */
     ComputeEquationResults computeStateEquation();
    /**
     * @brief Get the Components object
     * Returns a tuple with either the readableName (if available) or id, handle to it and whether it is a proxy
     * @return std::vector<std::tuple<str::string,RCPLIB::RCP<BGElement> ,bool>>& 
     */
     std::vector<std::tuple<std::string,const RCPLIB::RCP<BGElement> &,bool>> getComponents();
    /**
     * @brief Get the Element By Name object
     * Returns element with the specified name, throws Exception if not found
     * @param name 
     * @return const RCPLIB::RCP<BGElement>& 
     */
     virtual const RCPLIB::RCP<BGElement> & getElementByName(std::string name);
    /**
     * @brief Get the list of source elements associated with this bondgraph, entries are valid only after computeStateEquation is called
     * 
     * @return std::vector<RCPLIB::RCP<BGElement> >& List of elements
     */
     virtual std::vector<RCPLIB::RCP<BGElement> > &getSources();
    /**
     * @brief Get the list of Bonds associated with this bondgraph, entries are valid only after computeStateEquation is called
     * 
     * @return const std::vector<RCPLIB::RCP<BondInterface>>& List of bonds
     */
     virtual const std::vector<RCPLIB::RCP<BondInterface>> &getBonds() const;

    //IO related methods
    //!Serialize the bond graph in JSON form.
     std::string serialize(bool inIndent = false) const;
    //!Serialize the bond graph in flattend one component CellML form.
     std::string serializeAsCellML(std::string modelName);
     virtual void read(std::string filename);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraphInterface> &p);
    friend void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraphInterface> &p);
    friend  RCPLIB::RCP<BondGraphInterface>  createBondGraph();
};



/**
 * @brief Create a Resistane instance
 * Library users are not expected to call this as element definition is incomplete and domain is not assigned.
 * Users should call one of the generated factory methods
 * @return Reference counted pointer to the element 
 */
RCPLIB::RCP<BGElement> __createResistance_do_not_call_outside_of_lib(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

/**
 * @brief Create a Capacitance instance
 * Library users are not expected to call this as element definition is incomplete and domain is not assigned.
 * Users should call one of the generated factory methods 
 * @return Reference counted pointer to the element 
 */
RCPLIB::RCP<BGElement> __createCapacitance_do_not_call_outside_of_lib(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

/**
 * @brief Create a Inductance instance
 * Library users are not expected to call this as element definition is incomplete and domain is not assigned.
 * Users should call one of the generated factory methods
 * @return Reference counted pointer to the element 
 */
RCPLIB::RCP<BGElement> __createInductance_do_not_call_outside_of_lib(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);


} // namespace BG

