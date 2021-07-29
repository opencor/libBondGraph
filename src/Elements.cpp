#include "Bond.h"
#include "Elements.h"
#include "ElementsImpl.h"
#include "Exceptions.h"
#include "Port.h"
#include <sstream>
#include <symengine/parser.h>
#include <symengine/subs.h>
#include <symengine/visitor.h>
#include <units.hpp>

namespace BG {

#define SETIDANDUPDATE(idvar, target) \
    dofID = idvar; \
    std::ostringstream ss; \
    for (int ix = getNumStates(); ix < target.size(); ix++) { \
        if (target[ix]->universalConstant) \
            return; \
        ss.str(""); \
        ss.clear(); \
        ss << target[ix]->prefix << "_" << inId; \
        target[ix]->name = ss.str(); \
    }

Resistance::Resistance(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;
    mElementType = eResistance;
    mComponentGroup = eR;
}

void Resistance::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 0)
        throw BGException("Attempt to connect Single port component (Resistance) to more than one port!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void Resistance::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

Capacitance::Capacitance(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eCapacitance;
    mComponentGroup = eS;
}

void Capacitance::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 0)
        throw BGException("Attempt to connect Single port component (Capacitance) to more than one port!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void Capacitance::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

Inductance::Inductance(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eInductance;
    mComponentGroup = eS;
}

void Inductance::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 0)
        throw BGException("Attempt to connect Single port component (Inductance) to more than one port!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void Inductance::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

PotentialSource::PotentialSource(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = ePotentialSource;
    mComponentGroup = eU;
    RCPLIB::RCP<Value> parameter = RCPLIB::rcp(new Value("u", 1.0));
    mParameter.push_back(parameter);
    parameter->units = units::unit_from_string("V");
    numStates = 0;
    setName("Se");
    //Dont use e for effort, as symengine parse considers it as exponential; see Parser::parse_identifier of parser.cpp
    constitutiveEq = {"e_0 -u"};
    constitutiveEqIndex = {-1};
}

void PotentialSource::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 0)
        throw BGException("Attempt to connect Single port component (PotentialSource) to more than one port!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void PotentialSource::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

FlowSource::FlowSource(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eFlowSource;
    mComponentGroup = eU;
    RCPLIB::RCP<Value> parameter = RCPLIB::rcp(new Value("i", 1.0));
    mParameter.push_back(parameter);
    parameter->units = units::unit_from_string("A");
    numStates = 0;
    setName("Sf");
    //Dont use e for effort, as symengine parse considers it as exponential; see Parser::parse_identifier of parser.cpp
    constitutiveEq = {"f_0 -i"};
    constitutiveEqIndex = {-1};
}

void FlowSource::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 0)
        throw BGException("Attempt to connect Single port component (FlowSource) to more than one port!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void FlowSource::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

Concentration::Concentration(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = bConcentration;
    mComponentGroup = eS;
    mDomain = "Biochemical";
    RCPLIB::RCP<Value> mValue = RCPLIB::rcp(new Value("a_", 0.0)); // Molar quantity
    mParameter.push_back(mValue);
    RCPLIB::RCP<Value> mGasConstant = RCPLIB::rcp(new Value("R", 8.3145)); // Gas Constant
    mParameter.push_back(mGasConstant);
    RCPLIB::RCP<Value> mTemperature = RCPLIB::rcp(new Value("T", 293.14)); // Temperature
    mParameter.push_back(mTemperature);
    RCPLIB::RCP<Value> mBiochemicalConst = RCPLIB::rcp(new Value("k", 1.0)); // Biochemical Constant; exp(mu_0/RT)/V
    mParameter.push_back(mBiochemicalConst);

    numStates = 1;

    // SI units kg mol -1
    {
        auto un = units::unit_from_string("kg/mol");
        auto mult = un.multiplier();
        auto baseU = un.base_units();
        mValue->units = units::precise_unit(baseU, mult);
    }
    //Gas Constant kg⋅m2·K−1⋅mol−1s−2
    {
        auto un = units::unit_from_string("J/(K.mol)");
        auto mult = un.multiplier();
        auto baseU = un.base_units();
        mGasConstant->units = units::precise_unit(baseU, mult);
        mGasConstant->universalConstant = true;
    }
    //Temperature
    mTemperature->units = units::unit_from_string("K");
    mTemperature->universalConstant = true;
    //Biochemical Const
    {
        auto un = units::unit_from_string("mol/kg");
        auto mult = un.multiplier();
        auto baseU = un.base_units();
        mBiochemicalConst->units = units::precise_unit(baseU, mult);
    }
    setName("Ce");
    constitutiveEq = {"e_0 - R*T*log(k*a_0)", "f_0 - dot_a_0"};
    constitutiveEqIndex = {-1, -1};
}

void Concentration::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 0)
        throw BGException("Attempt to connect Single port component (Concentration) to more than one port!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void Concentration::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

Transformer::Transformer(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eTransformer;
    mComponentGroup = eJ;
    RCPLIB::RCP<Value> parameter = RCPLIB::rcp(new Value("n", 1.0));
    parameter->units = units::precise::one; //same as dimensionless
    mParameter.push_back(parameter);
    numStates = 0;

    setName("TF");
    constitutiveEq = {"e_1 - n * e_0", "f_0 + n * f_1"};
    constitutiveEqIndex = {-1, -1};
}

void Transformer::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 1)
        throw BGException("Attempt to connect Two port component (Transformer)  to more than two ports!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void Transformer::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

Gyrator::Gyrator(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eGyrator;
    mComponentGroup = eJ;
    RCPLIB::RCP<Value> parameter = RCPLIB::rcp(new Value("n", 1.0));
    parameter->units = units::precise::one; //same as dimensionless
    mParameter.push_back(parameter);
    numStates = 0;

    setName("GY");
    constitutiveEq = {"e_1 + n*f_0", "e_0 - n*f_1"};
    constitutiveEqIndex = {-1, -1};
}

void Gyrator::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > 1)
        throw BGException("Attempt to connect Two port component (Gyrator)  to more than two ports!");
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void Gyrator::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

ZeroJunction::ZeroJunction(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eZero;
    mComponentGroup = eJ;
    numStates = 0;
    portModified = true;
    setName("ZJ");
};

void ZeroJunction::connect(const RCPLIB::RCP<Port> &inPort)
{
    portModified = true;
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void ZeroJunction::disconnect(const RCPLIB::RCP<Port> &inPort)
{
    portModified = true;
    std::vector<RCPLIB::RCP<Port>>::iterator pHandle = std::find(mPorts.begin(), mPorts.end(), inPort);
    if (pHandle != mPorts.end())
        mPorts.erase(pHandle);
}

std::vector<std::string> &ZeroJunction::getConstitutiveEquations()
{
    if (portModified) {
        int numPorts = mPorts.size();
        std::ostringstream ss, fs;
        constitutiveEq.clear();
        constitutiveEqIndex.clear();
        ss << "e_" << numPorts - 1;
        std::string eEnd = ss.str();
        for (int i = numPorts - 2; i > -1; i--) {
            ss.str("");
            ss.clear();
            ss << "e_" << i << " - " << eEnd;
            constitutiveEq.push_back(ss.str());
            constitutiveEqIndex.push_back(-1);
            fs << "f_" << i << " + ";
        }
        fs << "f_" << numPorts - 1;

        constitutiveEq.push_back(fs.str());
        constitutiveEqIndex.push_back(-1);
        portModified = false;
    }
    return constitutiveEq;
}

/*
void ZeroJunction::setDof(size_t inId)
{
    mId = inId;
}
*/
OneJunction::OneJunction(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = eOne;
    mComponentGroup = eJ;
    numStates = 0;
    portModified = true;
    setName("OJ");
};

void OneJunction::connect(const RCPLIB::RCP<Port> &inPort)
{
    portModified = true;
    mPorts.push_back(inPort);
    inPort->connect(rcpPtr);
}

void OneJunction::disconnect(const RCPLIB::RCP<Port> &inPort)
{
    portModified = true;
    std::vector<RCPLIB::RCP<Port>>::iterator pHandle = std::find(mPorts.begin(), mPorts.end(), inPort);
    if (pHandle != mPorts.end())
        mPorts.erase(pHandle);
}

std::vector<std::string> &OneJunction::getConstitutiveEquations()
{
    if (portModified) {
        std::ostringstream ss, fs;
        constitutiveEq.clear();
        constitutiveEqIndex.clear();
        int numPorts = mPorts.size();
        double dsigma0 = mPorts[numPorts - 1]->getWeight();
        std::string sigma0 = "";
        if (dsigma0 != 1.0) {
            ss << dsigma0 << "*";
            sigma0 = ss.str();
        }
        ss.str("");
        ss.clear();
        ss << "f_" << numPorts - 1;
        std::string fEnd = ss.str();
        fs << sigma0 << "e_" << numPorts - 1;
        for (int i = numPorts - 2; i > -1; i--) {
            ss.str("");
            ss.clear();
            double sigma = mPorts[i]->getWeight();
            if (sigma != 1.0) {
                ss << sigma << "*";
                fs << " + " << sigma << "*e_" << i;
            } else {
                fs << " + e_" << i;
            }
            ss << "f_" << i << " - " << sigma0 << fEnd;
            auto p = SymEngine::parse(ss.str());
            ss.str("");
            ss.clear();
            ss << *p;
            constitutiveEq.push_back(ss.str());
            constitutiveEqIndex.push_back(-1);
        }
        auto p = SymEngine::parse(fs.str());
        fs.str("");
        fs.clear();
        fs << *p;
        constitutiveEq.push_back(fs.str());
        constitutiveEqIndex.push_back(-1);
        portModified = false;
    }
    return constitutiveEq;

    //Compute based on current connections, load and then return
    /*
            relations = []

            var = list(self._port_vectors().items())
            (e_0, f_0), port = var.pop()

            sigma_0 = port.weight
            partial_sum = sigma_0 * e_0

            while var:
                (e_i, f_i), port = var.pop()
                sigma_i = port.weight
                partial_sum += sigma_i * e_i
                relations.append(sigma_i * f_i - sigma_0 * f_0)

            relations.append(partial_sum)
            return relations

            */
}

/*
void OneJunction::setDof(size_t inId)
{
    mId = inId;
}
*/

Reaction::Reaction(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    mElementType = bReaction;
    mComponentGroup = eR;
    mDomain = "Biochemical";
    mParameter.clear();
    RCPLIB::RCP<Value> parameter = RCPLIB::rcp(new Value("r", 1.0)); //Reaction Rate
    mParameter.push_back(parameter);
    RCPLIB::RCP<Value> mGasConstant = RCPLIB::rcp(new Value("R", 8.3145)); // Gas Constant
    mParameter.push_back(mGasConstant);
    RCPLIB::RCP<Value> mTemperature = RCPLIB::rcp(new Value("T", 293.14)); // Temperature
    mParameter.push_back(mTemperature);
    numStates = 0;

    //Gas Constant kg⋅m2·K−1⋅mol−1s−2
    {
        auto un = units::unit_from_string("J/(K.mol)");
        auto mult = un.multiplier();
        auto baseU = un.base_units();
        mGasConstant->units = units::precise_unit(baseU, mult);
        mGasConstant->universalConstant = true;
    }

    //Temperature
    mTemperature->units = units::unit_from_string("K");
    mTemperature->universalConstant = true;
    //Reaction rate, First order
    {
        auto un = units::unit_from_string("mol/(s.L)");
        auto mult = un.multiplier();
        auto baseU = un.base_units();
        parameter->units = units::precise_unit(baseU, mult);
    }

    setName("Rx");
    constitutiveEq = {"f_0 + f_1", "f_0 - r*(exp(e_0/R/T) - exp(e_1/R/T))"};
    constitutiveEqIndex = {-1, -1};
    allocatedPorts = true;
    mPorts.clear();
    const RCPLIB::RCP<Port> lFromPort = createPort(false);
    const RCPLIB::RCP<Port> lToPort = createPort(true);
    mPorts.push_back(lFromPort);
    mPorts.push_back(lToPort);

}

void Reaction::connect(const RCPLIB::RCP<Port> &inPort)
{
    //inPort->connect(RCPLIB::rcp(this,true));
    inPort->connect(rcpPtr);
}

void Reaction::setDof(size_t inId) {
    SETIDANDUPDATE(inId, mParameter)}

//Userdefined bondgraph element, should fall within on of the major passivetype groups
UserDefined::UserDefined(const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;

    numStates = 0;
    maxPorts = 1;
    mElementType = eUserDefined;
    mComponentGroup = eU;
    RCPLIB::RCP<Value> parameter = RCPLIB::rcp(new Value("ue", "0.0"));
    mParameter.push_back(parameter);
    parameter->units = units::unit_from_string("V");
    RCPLIB::RCP<Value> parameter1 = RCPLIB::rcp(new Value("uf", "0.0"));
    mParameter.push_back(parameter1);
    parameter1->units = units::unit_from_string("A");

    numStates = 0;
    setName("Ux");
    //Dont use e for effort, as symengine parse considers it as exponential; see Parser::parse_identifier of parser.cpp
    constitutiveEq = {"e_0 - ue", "f_0 - uf"};
    constitutiveEqIndex = {-1, -1};
}

void UserDefined::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > maxPorts) {
        std::ostringstream ss;
        ss << "Userdefined " << maxPorts << "'s component connected to more than " << maxPorts << "!" << std::endl;
        throw BGException(ss.str());
    }
    if (maxPorts == 1) {
        mPorts.push_back(inPort);
        inPort->connect(rcpPtr);
    } else {
        inPort->connect(rcpPtr);
    }
}

void UserDefined::setDof(size_t inId)
{
    SETIDANDUPDATE(inId, mParameter)
}

void UserDefined::setPhysicalDomain(PhysicalDomain domain)
{
    if (portDimensions.find(domain) != portDimensions.end()) {
        auto res = portDimensions[domain];
        mParameter[0]->units = units::unit_from_string(std::get<0>(res));
        mParameter[1]->units = units::unit_from_string(std::get<1>(res));
    } else {
        throw BGException("Unknown physical domain " + domain);
    }
}

void UserDefined::setConstitutiveEqIndex(int index, int eqNo)
{
    //throw BGException("Cannot set constitutive equations for UserDefined types");
}

PortHamiltonian::PortHamiltonian(std::string hexpr,
                                 std::vector<std::string> states,
                                 const RCPLIB::RCP<ElementImpl>& data, std::string inId, bool proxy_)
    : BondGraphElementBase(data, inId)
{
    proxy = proxy_;
    if (proxy)
        return;
    numStates = 0;
    maxPorts = 0;
    mElementType = ePortHamiltonian;
    mComponentGroup = ePH;
    setName("PH");
    //Parse to get the states, rates and the parameters
    SymEngine::RCP<const SymEngine::Basic> expr = SymEngine::parse(hexpr);
    auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*expr);
    std::unordered_map<std::string, bool> derivatives;
    std::unordered_map<std::string, std::tuple<std::string, std::string>> stateNameMaps;
    std::vector<std::string> parameters;
    SymEngine::map_basic_basic subs;
    //All state names are suffixed with "_" for allowing global dof naming
    for (int i = 0; i < states.size(); i++) {
        std::string st = states[i];
        std::string res = st + "_";
        derivatives[res] = false;
        subs[SymEngine::parse(st)] = SymEngine::parse(st + "_");
    }

    for (auto at : atoms) {
        std::string sym = at->__str__();
        if (std::find(states.begin(), states.end(), sym) == states.end()) {
            if (sym.find("dot_", 0) == 0) {
                std::string statevar = sym.substr(4) + "_";
                if (derivatives.find(statevar) != derivatives.end()) {
                    derivatives[statevar] = true;
                    subs[SymEngine::parse(sym)] = SymEngine::parse(sym + "_");
                } else
                    logWarn("Derivative ",sym," does not have a prescribed state for port hamiltonian ",hexpr);
            } else {
                parameters.push_back(sym);
            }
        }
    }
    auto nexpr = expr->subs(subs);
    logInfo("Transformed port hamiltonian expression ",expr->__str__(),"to",nexpr->__str__());
    std::ostringstream ss;
    for (int i = 0; i < states.size(); i++) {
        std::string st = states[i];
        std::string coord = st + "_";
        ss.str("");
        ss.clear();

        //Potential = rate of change of Hamiltonian(Energy) =  \partial{H}/\partial{st}
        //Flow = rate of change of coordinate
        auto e_0 = nexpr->diff(SymEngine::symbol(coord));
        ss<<"e_"<<i<<"-"<<*e_0;
        constitutiveEq.push_back(ss.str());
        constitutiveEqIndex.push_back(-1);
        constitutiveEq.push_back("dot_"+st+" - f_"+std::to_string(i));
        constitutiveEqIndex.push_back(-1);
        numStates++;
        RCPLIB::RCP<Value> mValue = RCPLIB::rcp(new Value(coord, 0.0));
        mValue->units = units::precise::one; //Same as dimensionless
        mParameter.push_back(mValue);
    }
    maxPorts = numStates; //Number of ports = number of states
    for (auto p : parameters) {
        RCPLIB::RCP<Value> mValue = RCPLIB::rcp(new Value(p, 0.0));
        mValue->units = units::precise::one; //Same as dimensionless
        mParameter.push_back(mValue);
    }
    //Determine which components need gyrators to ensure effort-flow consistency, for instamce if the
    //Hamiltonion has a both electrical and mechanical domains, then for the electrical storage effort and flow are as expected
    //However, for the mechnical one they are inverted
    allocatedPorts = true;
    mPorts.clear();
    for(int i=0;i<maxPorts;i++){
        mPorts.push_back(createPort(false));
    }
}

void PortHamiltonian::connect(const RCPLIB::RCP<Port> &inPort)
{
    if (mPorts.size() > maxPorts) {
        std::ostringstream ss;
        ss << "Userdefined " << maxPorts << "'s component connected to more than " << maxPorts << "!" << std::endl;
        throw BGException(ss.str());
    }
    if (maxPorts == 1) {
        mPorts.push_back(inPort);
        inPort->connect(rcpPtr);
    } else {
        inPort->connect(rcpPtr);
    }
}

void PortHamiltonian::setDof(size_t inId)
{
    SETIDANDUPDATE(inId, mParameter)
}

void PortHamiltonian::setPhysicalDomain(PhysicalDomain domain)
{
    throw BGException("Port hamiltonian can span multiple physical domains, use symmetric gyrator(s) to ensure flow/effort coupling consistency ");
}

void PortHamiltonian::setConstitutiveEqIndex(int index, int eqNo)
{
    constitutiveEqIndex[eqNo] = index;
    //throw BGException("Cannot set constitutive equations for Port Hamiltonian types");
}

void PortHamiltonian::addConstraints(std::string cons)
{
    auto trim = [](const std::string s) {
        auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
        auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
        return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
    };

    if (std::find(constraints.begin(), constraints.end(), trim(cons)) == constraints.end()) {
        constraints.push_back(trim(cons));
    }
}

std::vector<std::string> &PortHamiltonian::getConstraints()
{
    return constraints;
}


std::string PortHamiltonian::getHamiltonian()
{
    return hamiltonian;
}

} // namespace BG
