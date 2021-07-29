#ifndef ELEMENTS_H_
#define ELEMENTS_H_
#include "Elementsbase.h"
#include "friends.h"
#include "Port.h"
#include "RCP.h"
#include "export.h"
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

//Forward declaration
class Port;
class Bond;

class EXPORTED Resistance: public BondGraphElementBase
{
    Resistance(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Resistance> __createResistance_do_not_call_outside_of_lib(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
    DEFINE_FRIENDS_OF_RESISTANCE
};

class EXPORTED Capacitance: public BondGraphElementBase
{
    Capacitance(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Capacitance> __createCapacitance_do_not_call_outside_of_lib(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
    DEFINE_FRIENDS_OF_CAPACITANCE
};

class EXPORTED Inductance: public BondGraphElementBase
{
    Inductance(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Inductance> __createInductance_do_not_call_outside_of_lib(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
    DEFINE_FRIENDS_OF_INDUCTANCE
};

class EXPORTED PotentialSource: public BondGraphElementBase
{
    PotentialSource(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<PotentialSource> createPotentialSource(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
    DEFINE_FRIENDS_OF_PotentialSource
};

class EXPORTED FlowSource: public BondGraphElementBase
{
    FlowSource(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<FlowSource> createFlowSource(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
    DEFINE_FRIENDS_OF_FLOWSOURCE
};

class EXPORTED Transformer: public BondGraphElementBase
{
    Transformer(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Transformer> createTransformer(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED Gyrator: public BondGraphElementBase
{
    Gyrator(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Gyrator> createGyrator(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED Concentration: public BondGraphElementBase
{
    Concentration(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Concentration> createConcentration(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED Reaction: public BondGraphElementBase
{
    Reaction(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    friend RCPLIB::RCP<Reaction> createReaction(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED PortHamiltonian: public BondGraphElementBase
{
    PortHamiltonian(std::string expr, 
                    std::vector<std::string> states,
                    const RCPLIB::RCP<ElementImpl>& data, 
                    std::string inId = "-1", bool proxy = false);
    unsigned int maxPorts;
    std::vector<std::string> constraints;
public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    
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
    friend RCPLIB::RCP<PortHamiltonian> createPortHamiltonian(std::string expr,
                                                              std::vector<std::string> states, const RCPLIB::RCP<ElementImpl>& data);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraphElementBase> &p);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED UserDefined: public BondGraphElementBase
{
    UserDefined(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);
    unsigned int maxPorts;

public:
    //! Connect the component to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    //virtual void setId(std::string inId);
    virtual void setDof(size_t inId);
    //! Set constitutive equation
    void setConstitutiveEqIndex(int index, int eqNo = 0);
    //! Set the physical domain, this will not change constitutive relations but set the port dimensions
    void setPhysicalDomain(PhysicalDomain domain);
    friend RCPLIB::RCP<UserDefined> createUserDefined(const RCPLIB::RCP<ElementImpl>& data);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraphElementBase> &p);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED ZeroJunction: public BondGraphElementBase
{
    ZeroJunction(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    virtual ~ZeroJunction()
    {
    }

    PassiveType getType() const
    {
        return mElementType;
    }
    //! Connect the junction to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    void disconnect(const RCPLIB::RCP<Port> &inPort);
    std::vector<std::string> &getConstitutiveEquations();
    //virtual void setId(std::string inId);
    friend RCPLIB::RCP<ZeroJunction> createZeroJunction(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};

class EXPORTED OneJunction: public BondGraphElementBase
{
    OneJunction(const RCPLIB::RCP<ElementImpl>& data, std::string inId = "-1", bool proxy = false);

public:
    virtual ~OneJunction()
    {
    }

    PassiveType getType() const
    {
        return mElementType;
    }
    //! Connect the junction to a port.
    void connect(const RCPLIB::RCP<Port> &inPort);
    void disconnect(const RCPLIB::RCP<Port> &inPort);
    std::vector<std::string> &getConstitutiveEquations();
    //virtual void setId(std::string inId);
    friend RCPLIB::RCP<OneJunction> createOneJunction(const RCPLIB::RCP<ElementImpl>& data);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &data);
};
} // namespace BG

#endif