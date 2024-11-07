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
#include "Elementsbase.h"
#include "friends.h"
#include <algorithm>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>


namespace BG {

/**
     * @brief Concrete implementations of standard bond graph elements
     * 
     */

class  Resistance: public BondGraphElementBase
{
    Resistance(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<BGElement> __createResistance_do_not_call_outside_of_lib(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
    DEFINE_FRIENDS_OF_RESISTANCE
};

class  Capacitance: public BondGraphElementBase
{
    Capacitance(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<BGElement> __createCapacitance_do_not_call_outside_of_lib(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
    DEFINE_FRIENDS_OF_CAPACITANCE
};

class  Inductance: public BondGraphElementBase
{
    Inductance(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<BGElement> __createInductance_do_not_call_outside_of_lib(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
    DEFINE_FRIENDS_OF_INDUCTANCE
};

class  PotentialSource: public BondGraphElementBase
{
    PotentialSource(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createPotentialSource(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
    DEFINE_FRIENDS_OF_EFFORTSOURCE
};

class  FlowSource: public BondGraphElementBase
{
    FlowSource(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createFlowSource(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
    DEFINE_FRIENDS_OF_FLOWSOURCE
};

class  Transformer: public BondGraphElementBase
{
    Transformer(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createTransformer(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  Gyrator: public BondGraphElementBase
{
    Gyrator(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createGyrator(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  Concentration: public BondGraphElementBase
{
    Concentration(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createConcentration(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  Reaction: public BondGraphElementBase
{
    Reaction(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createReaction(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  ChemicalSource: public BondGraphElementBase
{
    ChemicalSource(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createChemostat(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  ChemicalFlux: public BondGraphElementBase
{
    ChemicalFlux(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createFlowstat(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};


class Stoichiometry: public BondGraphElementBase
{
    Stoichiometry(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    virtual void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend  RCPLIB::RCP<BGElement> createStoichiometry(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  PortHamiltonian: public BondGraphElementBase
{
    PortHamiltonian(std::string expr, 
                    std::vector<std::string> states,
                    const RCPLIB::RCP<BGElementData>& data, 
                    std::string inId = "-1", bool proxy = false);
    unsigned int maxPorts;
    std::vector<std::string> constraints;
public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    
    virtual void setDof(size_t inId);
     std::string getHamiltonian();
    //! In case of port hamiltonians these are interpreted as constraints
    void setConstitutiveEqIndex(int index, int eqNo = 0);
    //! Add constraints
     void addConstraints(std::string cons);
    //! Get list of constraints
     std::vector<std::string>& getConstraints();
    //! Set the physical domain, this will not change constitutive relations but set the port dimensions
     void setPhysicalDomain(PhysicalDomain domain);
    friend  RCPLIB::RCP<BGElement> createPortHamiltonian(std::string expr,
                                                              std::vector<std::string> states, const RCPLIB::RCP<BGElementData>& data);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BGElement>  &p);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  UserDefined: public BondGraphElementBase
{
    UserDefined(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);
    unsigned int maxPorts;

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    //! Set constitutive equation
    void setConstitutiveEqIndex(int index, int eqNo = 0);
    //! Set the physical domain, this will not change constitutive relations but set the port dimensions
     void setPhysicalDomain(PhysicalDomain domain);
    friend  RCPLIB::RCP<BGElement> createUserDefined(const RCPLIB::RCP<BGElementData>& data);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BGElement>  &p);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  ZeroJunction: public BondGraphElementBase
{
    ZeroJunction(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    virtual ~ZeroJunction()
    {
    }

    PassiveType getType() const
    {
        return mElementType;
    }
    //! Connect the junction to a port.
    virtual  void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    virtual  void disconnect(const RCPLIB::RCP<PortInterface> &inPort);
     std::vector<std::string> &getConstitutiveEquations();
    //virtual void setId(std::string inId);
    friend  RCPLIB::RCP<BGElement> createZeroJunction(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};

class  OneJunction: public BondGraphElementBase
{
    OneJunction(const RCPLIB::RCP<BGElementData>& data, std::string inId = "-1", bool proxy = false);

public:
    virtual ~OneJunction()
    {
    }

    PassiveType getType() const
    {
        return mElementType;
    }
    //! Connect the junction to a port.
    virtual  void connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum=0);
    virtual  void disconnect(const RCPLIB::RCP<PortInterface> &inPort);
     std::vector<std::string> &getConstitutiveEquations();
    //virtual void setId(std::string inId);
    friend  RCPLIB::RCP<BGElement> createOneJunction(const RCPLIB::RCP<BGElementData>& data);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &data);
};
} // namespace BG
