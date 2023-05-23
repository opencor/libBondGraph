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
#include "bondgraph.hpp"
#include "friends.h"
#include "export.h"
#include <string>
#include <units.hpp>
#include <unordered_map>
#include <vector>

namespace BG {


/*! \brief Bond Graph BondGraphElement
    *  Encapsulates a bond graph element. 
    */

class  BondGraphElementBase: public BGElement
{
protected:
    RCPLIB::RCP<BGElementData> _data;
    std::string parent = "";
    bool proxy = false; //! True if this is a proxy, BGElementData data is not owned by this instance, the owner has the same id

    std::string &mName; //!< Name of the element
    std::string &hamiltonian; //!< Hamltonian energy description (valid if the element is a Port Hamiltonian)
    std::string &mId; //!< Element's ID
    PassiveType &mElementType; //!< Type of element
    ComponentGroup &mComponentGroup; //!< Component group association.
    PhysicalDomain &mDomain; //! Physical Domain
    unsigned int &numStates; //Number of states for this Bondgraph element
    bool &allocatedPorts; //! True if ports are preallocated
    bool &portModified; //! True if port assignment was modified, used by junctions
    std::vector<RCPLIB::RCP<PortInterface>> &mPorts; //!< Ports connected to the component. Size depends on the type of component.
    std::vector<std::string> &constitutiveEq; //!< Constitutive equations
    std::vector<int> &constitutiveEqIndex; //!< Index of constitutive equation in connectivity matrix
    std::vector<RCPLIB::RCP<Value>> &mParameter; //!Array of state and parameter values associated with this bondgraph element
    RCPLIB::RCP<BGElement>  rcpPtr; //RCP pointer for this
    size_t dofID = 0; //! Dof id  for the element
    //Map of physical domains
    DEFINE_PHYSICAL_DOMAINS
    BondGraphElementBase(RCPLIB::RCP<BGElementData> data, std::string inId = "-1");

    //! Proxy constructor, the data is shared and not created
    BondGraphElementBase(const RCPLIB::RCP<BGElement>  &data);

public:
    virtual ~BondGraphElementBase();
    //! Get element name.
     const std::string &getName() const;
    //! Set element name
     void setName(const std::string &name);
    //! Get variable name for this element
    std::string getVariableName();
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
     std::vector<RCPLIB::RCP<PortInterface>> &getPorts();
    //! Return port at pin i
     RCPLIB::RCP<PortInterface> &getPorts(unsigned int i);
    //! Connect the component of a port
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0) = 0;
    //! Connect the compomemy to inPort and disconnect it from outPort
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, const RCPLIB::RCP<PortInterface> &outPort);
    //! Disconnect the component of a port
     virtual void disconnect(const RCPLIB::RCP<PortInterface> &inPort);
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
    //! Update parameter name
    virtual bool updateParameterName(std::string current,std::string target);
    friend  RCPLIB::RCP<BGElement>  createProxy(const RCPLIB::RCP<BGElement>  &obj);
    friend  RCPLIB::RCP<BGElement>  createClone(const RCPLIB::RCP<BGElement>  &obj);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BGElement>  &p);
    friend void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraphInterface> &p);
};

template<typename T>
inline RCPLIB::RCP<Value> BondGraphElementBase::setValue(std::string name, T value)
{
    return setValue(name,std::to_string(value));
}

template<typename T>
inline RCPLIB::RCP<Value> BondGraphElementBase::setParameter(std::string name, T value, std::string unit)
{
    return setParameter(name, std::to_string(value),unit);
}

std::string replaceAll(std::string str, const std::string &from, const std::string &to);
} // namespace BG

