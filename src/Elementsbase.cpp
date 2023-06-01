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


#include "Bond.h"
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

BondGraphElementBase::BondGraphElementBase(RCPLIB::RCP<BGElementData> data, std::string inId)
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

BondGraphElementBase::BondGraphElementBase(const RCPLIB::RCP<BGElement>  &data_)
    : _data(RCPLIB::null),
    mId(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mId),
    mName(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mName),
    hamiltonian(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->hamiltonian),
    mElementType(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mElementType),
    mComponentGroup(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mComponentGroup),
    mDomain(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mDomain),
    numStates(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->numStates),
    allocatedPorts(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->allocatedPorts),
    portModified(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->portModified),

    mPorts(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mPorts),
    constitutiveEq(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->constitutiveEq),
    constitutiveEqIndex(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->constitutiveEqIndex),
    mParameter(RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->mParameter)
{
    parent = RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(data_)->parent;
    annotation = {};
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
//!Get variable name
std::string BondGraphElementBase::getVariableName(){
    std::string varname = mName;
    if(isdigit(mName[0])){
        if(mName[0]=='0'){
        varname[0] = 'O';
        }else if(mName[0]=='1'){
        varname[0] = 'I';
        }else{
        varname = "E_"+mName;
        }
    }
    //Handle :
    std::replace(varname.begin(), varname.end(), ':', 'c');
    std::replace(varname.begin(), varname.end(), '^', 'p');
    std::replace(varname.begin(), varname.end(), '{', 'o');
    std::replace(varname.begin(), varname.end(), '}', 'x');    
    return varname;
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
std::vector<RCPLIB::RCP<PortInterface>> &BondGraphElementBase::getPorts()
{
    return mPorts;
}
RCPLIB::RCP<PortInterface> &BondGraphElementBase::getPorts(unsigned int i)
{
    return mPorts[i];
}

// void BondGraphElementBase::connect(const RCPLIB::RCP<PortInterface> &inPort, int portNum)
// {
//     throw BGException("Connect not implemented at BondGraphElementBase level");
// }

void BondGraphElementBase::connect(const RCPLIB::RCP<PortInterface> &inPort, const RCPLIB::RCP<PortInterface> &outPort)
{
    std::vector<RCPLIB::RCP<PortInterface>>::iterator pHandle = std::find(mPorts.begin(), mPorts.end(), outPort);
    if (pHandle != mPorts.end()) {
        pHandle->release();
        //mPorts.erase(pHandle);
        mPorts[pHandle - mPorts.begin()] = inPort;
    }
    //mPorts.push_back(inPort);
    portModified = true;
}

void BondGraphElementBase::disconnect(const RCPLIB::RCP<PortInterface> &inPort)
{
    std::vector<RCPLIB::RCP<PortInterface>>::iterator pHandle = std::find(mPorts.begin(), mPorts.end(), inPort);
    if (pHandle != mPorts.end()) {
        pHandle->release();
        //mPorts.erase(pHandle);
    }
    portModified = true;
}

bool BondGraphElementBase::updateParameterName(std::string current,std::string target){
    for(auto& v : mParameter){
        if(v->name==current){
           for(int ci=0;ci<constitutiveEq.size();ci++){
               auto ce = constitutiveEq[ci];
               constitutiveEq[ci] = replaceAll(ce,v->name,target);
            }
            v->name = target;
            return true;
        }
    }
    return false;
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




RCPLIB::RCP<Value> BondGraphElementBase::setParameter(std::string name, std::string value, std::string unit_)
{
    bool found = false;
    std::string unit = unit_;
    auto un = units::unit_from_string(unit);
    auto mult = un.multiplier();
    auto baseU = un.base_units();
    auto preciseunit = units::to_string(units::precise_unit(baseU, mult));
    RCPLIB::RCP<Value> parameter;
    for (auto &param : mParameter) {
        if (name.rfind(param->prefix, 0) == 0) {
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
        auto preciseunit = units::to_string(units::precise::one); //Same as dimensionless or use units::unit_from_string("dimless");
        if(unit!=""){
            auto un = units::unit_from_string(unit);
            auto mult = un.multiplier();
            auto baseU = un.base_units();
            preciseunit = units::to_string(units::precise_unit(baseU, mult));
        }
        for (auto &param : mParameter) {
            if (name.rfind(param->prefix,0)==0) {
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



std::string BondGraphElementBase::getHamiltonian() {
    throw BGException("Not implemented for basic bondgraph types");
}

RCPLIB::RCP<Value> BondGraphElementBase::setValue(std::string name, std::string value)
{
    for (auto &param : mParameter) {
        if (name.rfind(param->prefix,0)==0) {
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

//! Set the PMR annotation
void BondGraphElementBase::setPMRAnnotation(nlohmann::json& annotation_){
    annotation = nlohmann::json(annotation_);
}

//! Get the PMR annotation
nlohmann::json& BondGraphElementBase::getPMRAnnotation(){
    return annotation;
}

} // namespace BG
