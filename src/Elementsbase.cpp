#include "Bond.h"
#include "ElementsImpl.h"
#include "Elementsbase.h"
#include "Exceptions.h"
#include "Port.h"
#include "logging.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <string>

namespace BG {

std::string replaceAll(std::string str, const std::string &from, const std::string &to)
{
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

std::ostream &operator<<(std::ostream &os, const Value &val)
{
    os << val.name << " " << val.value;
    return os;
}

BondGraphElementBase::BondGraphElementBase(RCPLIB::RCP<ElementImpl> data, std::string inId)
    : _data(data)
    , mId(data->mId)
    , mName(data->mName)
    , hamiltonian(data->hamiltonian)
    , mElementType(data->mElementType)
    , mComponentGroup(data->mComponentGroup)
    , mDomain(data->mDomain)
    , numStates(data->numStates)
    , allocatedPorts(data->allocatedPorts)
    , portModified(data->portModified)
    , mPorts(data->mPorts)
    , constitutiveEq(data->constitutiveEq)
    , constitutiveEqIndex(data->constitutiveEqIndex)
    , mParameter(data->mParameter)
{
    std::ostringstream ss;
    if (inId != "-1") {
        ss << inId;
    } else {
        ss << rand();
    }
    mId = ss.str();
    logDebug("Created element ", mName, " with id ", mId, " type ", mElementType);
}

BondGraphElementBase::BondGraphElementBase(const RCPLIB::RCP<BondGraphElementBase> &data)
    : _data(RCPLIB::null)
    , mId(data->mId)
    , mName(data->mName)
    , hamiltonian(data->hamiltonian)
    , mElementType(data->mElementType)
    , mComponentGroup(data->mComponentGroup)
    , mDomain(data->mDomain)
    , numStates(data->numStates)
    , allocatedPorts(data->allocatedPorts)
    , portModified(data->portModified)

    , mPorts(data->mPorts)
    , constitutiveEq(data->constitutiveEq)
    , constitutiveEqIndex(data->constitutiveEqIndex)
    , mParameter(data->mParameter)
{
    parent = data->parent;
}

BondGraphElementBase::~BondGraphElementBase()
{

}

//! Get element name.
const std::string &BondGraphElementBase::getName() const
{
    return mName;
}
//! Set element name
void BondGraphElementBase::setName(const std::string &name)
{
    mName = name;
}

//! Return the element ID.
std::string BondGraphElementBase::getId() const
{
    return mId;
}
//! Set the element ID.
void BondGraphElementBase::setId(std::string inId)
{
    mId = inId;
}

//! Set the DoF ID.
void BondGraphElementBase::setDof(size_t inId)
{
    dofID = inId;
}
//! Get the Dof ID.
size_t BondGraphElementBase::getDof()
{
    return dofID;
}

void BondGraphElementBase::setDomain(PhysicalDomain domain)
{
    mDomain = domain;
}

PhysicalDomain BondGraphElementBase::getDomain()
{
    return mDomain;
}

//! Set the element ID.
void BondGraphElementBase::setId(long int inId)
{
    mId = std::to_string(inId);
}
//! Set the parent BondGraph's name
void BondGraphElementBase::setParent(std::string &pid)
{
    parent = pid;
}
//! Get the parent BondGraph's name
const std::string BondGraphElementBase::getParent()
{
    return parent;
}

void BondGraphElementBase::setProxy(bool flag)
{
    proxy = flag;
}

bool BondGraphElementBase::isProxy()
{
    return proxy;
}

//Get the number of states for this bond graph element
unsigned int BondGraphElementBase::getNumStates()
{
    return numStates;
}

bool BondGraphElementBase::preallocatedPorts()
{
    return allocatedPorts;
}

void BondGraphElementBase::setPortModified(bool flag)
{
    portModified = flag;
}

//! Return the ports associated to the component.
std::vector<RCPLIB::RCP<Port>> &BondGraphElementBase::getPorts()
{
    return mPorts;
}
RCPLIB::RCP<Port> &BondGraphElementBase::getPorts(unsigned int i)
{
    return mPorts[i];
}

void BondGraphElementBase::connect(const RCPLIB::RCP<Port> &inPort)
{
    throw BGException("Connect not implemented at BondGraphElementBase level");
}

void BondGraphElementBase::connect(const RCPLIB::RCP<Port> &inPort, const RCPLIB::RCP<Port> &outPort)
{
    std::vector<RCPLIB::RCP<Port>>::iterator pHandle = std::find(mPorts.begin(), mPorts.end(), outPort);
    if (pHandle != mPorts.end()) {
        pHandle->release();
        mPorts.erase(pHandle);
    }
    mPorts.push_back(inPort);
    portModified = true;
}

void BondGraphElementBase::disconnect(const RCPLIB::RCP<Port> &inPort)
{
    std::vector<RCPLIB::RCP<Port>>::iterator pHandle = std::find(mPorts.begin(), mPorts.end(), inPort);
    if (pHandle != mPorts.end()) {
        pHandle->release();
        mPorts.erase(pHandle);
    }
    portModified = true;
}

std::vector<std::tuple<std::string, RCPLIB::RCP<Value>>> BondGraphElementBase::values()
{
    std::vector<std::tuple<std::string, RCPLIB::RCP<Value>>> mp;
    std::ostringstream ss;
    for (unsigned int i = 0; i < numStates; i++) {
        ss.str("");
        ss.clear();
        ss << mParameter[i]->prefix << i; //use prefix as name will change with each computeEquation call, prefix is local name
        mp.push_back(std::make_tuple(ss.str(), mParameter[i]));
    }
    for (unsigned int i = numStates; i < mParameter.size(); i++) {
        mp.push_back(std::make_tuple(mParameter[i]->prefix, mParameter[i]));
    }
    return mp;
};

template<typename T>
RCPLIB::RCP<Value> BondGraphElementBase::setParameter(std::string name, T value, std::string unit)
{
    return setParameter(name, std::to_string(value),unit);
}


RCPLIB::RCP<Value> BondGraphElementBase::setParameter(std::string name, std::string value, std::string unit)
{
    bool found = false;
    auto un = units::unit_from_string(unit);
    auto mult = un.multiplier();
    auto baseU = un.base_units();
    auto preciseunit = units::precise_unit(baseU, mult);
    RCPLIB::RCP<Value> parameter;
    for (auto &param : mParameter) {
        if (param->prefix == name) {
            parameter = param;
            param->value = value;
            param->units = preciseunit;
            found = true;
            break;
        }
    }
    if (!found) {
        parameter = RCPLIB::rcp(new Value(name, value));
        parameter->units = preciseunit;
        mParameter.push_back(parameter);
    }
    return parameter;
}

void BondGraphElementBase::setSIUnit(std::string name, std::string unit){
    if (mComponentGroup == ePH || mComponentGroup == eJ) {
        bool found = false;
        auto preciseunit = units::precise::one; //Same as dimensionless or use units::unit_from_string("dimless");
        if(unit!=""){
            auto un = units::unit_from_string(unit);
            auto mult = un.multiplier();
            auto baseU = un.base_units();
            preciseunit = units::precise_unit(baseU, mult);
        }
        for (auto &param : mParameter) {
            if (param->prefix == name) {
                param->units = preciseunit;
                found = true;
                break;
            }
        }
        if(!found){
            throw BGException("State/Parameter with name "+name+" not found!");
        }
    } else {
        throw BGException("SI Units can be set only for non Port Hamiltonian, transformer, gyrator elements !");
    }
}

template<typename T>
RCPLIB::RCP<Value> BondGraphElementBase::setValue(std::string name, T value)
{
    return setValue(name,std::to_string(value));
}

std::string BondGraphElementBase::getHamiltonian() {
    throw BGException("Not implemented for basic bondgraph types");
}

RCPLIB::RCP<Value> BondGraphElementBase::setValue(std::string name, std::string value)
{
    for (auto &param : mParameter) {
        if (param->prefix == name) {
            param->value = value;
            return param;
        }
    }
    throw BGException(name + " parameter/state not found.");
}

PassiveType BondGraphElementBase::getType() const
{
    return mElementType;
}
//! Return the component associated group.
ComponentGroup BondGraphElementBase::getComponentGroup() const
{
    return mComponentGroup;
}
//! Set constitutive equations index
void BondGraphElementBase::setConstitutiveEqIndex(int index, int eqNo)
{
    constitutiveEqIndex[eqNo] = index;
}
//! Get constitutive equations index
int BondGraphElementBase::getConstitutiveEqIndex(int eqNo)
{
    return constitutiveEqIndex[eqNo];
}
//! Get constitutive equations
std::vector<std::string> &BondGraphElementBase::getConstitutiveEquations()
{
    return constitutiveEq;
}

RCPLIB::RCP<Value> BondGraphElementBase::setUniversalConstant(std::string name, double &value, std::string unit)
{
    auto uCont = setParameter(name, value, unit);
    uCont->universalConstant = true;
    return uCont;
}

} // namespace BG
