#include "Bond.h"
#include "ElementsImpl.h"
#include "Exceptions.h"
#include "Port.h"
#include "bondgraph.h"
#include "componentregistry.h"
#include "thirdparty/json.hpp"
#include <units.hpp>

namespace BG {
//Forward declaration
RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &proxy = RCPLIB::null);

void to_json(nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &p)
{
    j["class"] = "ElementImpl";
    j["mName"] = p->mName;
    j["mId"] = p->mId;
    j["hamiltonian"] = p->hamiltonian;
    j["mElementType"] = p->mElementType;
    j["mComponentGroup"] = p->mComponentGroup;
    j["mDomain"] = p->mDomain;
    j["numStates"] = p->numStates;
    j["allocatedPorts"] = p->allocatedPorts;
    j["portModified"] = p->portModified;

    nlohmann::json ports;
    for (int i = 0; i < p->mPorts.size(); i++) {
        std::string ix = std::to_string(i);
        ports[ix] = p->mPorts[i];
    }
    j["mPorts"] = ports;
    j["constitutiveEq"] = p->constitutiveEq;
    j["constitutiveEqIndex"] = p->constitutiveEqIndex;
    nlohmann::json pms;
    for (int i = 0; i < p->mParameter.size(); i++) {
        std::string ix = std::to_string(i);
        pms[ix] = p->mParameter[i];
    }
    j["mParameter"] = pms;
}

void to_json(nlohmann::json &j, const RCPLIB::RCP<Port> &p)
{
    j["class"] = "Port";
    if (!p->mComponent.is_null()) {
        j["mComponent"] = p->mComponent->getId(); //!< Attached component
        auto &ports = p->mComponent->getPorts();
        j["portNo"] = std::find(ports.begin(), ports.end(), p) - ports.begin();
    } else {
        j["mComponent"] = NULL;
        j["portNo"] = NULL;
    }
    if (!p->mBond.is_null())
        j["mBond"] = p->mBond->mId; //!< Attached bond
    else
        j["mBond"] = NULL;
    j["mdofIndex"] = p->mdofIndex; //!< Index into dof dimensions
    j["weight"] = p->weight; //!< Port weight -1.0 for inverting
    j["receivesPower"] = p->receivesPower; //! True if power arrow is incident on this port
}

void to_json(nlohmann::json &j, const RCPLIB::RCP<Value> &p)
{
    j["class"] = "Value";
    j["prefix"] = p->prefix;
    j["name"] = p->name;
    j["value"] = p->value;
    j["units"] = units::to_string(p->units);
    j["universalConstant"] = p->universalConstant;
}

void to_json(nlohmann::json &j, const RCPLIB::RCP<Bond> &p)
{
    j["class"] = "Bond";
    j["mId"] = p->mId;
    j["mName"] = p->mName;
    j["mFromPort"] = p->mFromPort; //!< Origin port of the bond.
    j["mToPort"] = p->mToPort; //!< Destination port of the bond.

    j["mStorageType"] = p->mStorageType; //!< Storage type to complete the group type.
    j["mGroup"] = p->mGroup; //!< Bond group association.

    j["mStorageIndex"] = p->mStorageIndex;
    j["bid"] = p->bid;
}

void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraphElementBase> &p)
{
    j["class"] = "BondGraphElementBase";
    j["data"] = p->_data;
    j["parent"] = p->parent;
    j["proxy"] = p->proxy;
    j["dofID"] = p->dofID;
    j["readableName"] = "";

    const std::unordered_map<std::string, std::string> &nameMap = ComponentRegistry::defaultRegistry().getNameMap();
    if (nameMap.find(p->mId) != nameMap.end())
        j["readableName"] = nameMap.at(p->mId);

    if (p->mElementType == ePortHamiltonian) {
        const RCPLIB::RCP<PortHamiltonian> &ph = RCPLIB::rcp_dynamic_cast<PortHamiltonian>(p);
        j["maxPorts"] = ph->maxPorts;
        j["constraints"] = ph->constraints;
    }
    if (p->mElementType == eUserDefined) {
        const RCPLIB::RCP<UserDefined> &ud = RCPLIB::rcp_dynamic_cast<UserDefined>(p);
        j["maxPorts"] = ud->maxPorts;
    }
}

void to_json(nlohmann::json &j, const RCPLIB::RCP<BondGraph> &p)
{
    //First do the child bond graphs
    j["class"] = "BondGraph";
    j["id"] = p->id;
    j["cloned"] = p->cloned;
    j["ownerID"] = p->ownerID;
    j["readableName"] = "";
    const std::unordered_map<std::string, std::string> &nameMap = ComponentRegistry::defaultRegistry().getNameMap();
    if (nameMap.find(p->id) != nameMap.end())
        j["readableName"] = nameMap.at(p->id);
    std::vector<nlohmann::json> partsDesc;
    std::map<int, std::string> partOrder;
    for (int i = 0; i < p->mParts.size(); i++) {
        partOrder[i] = p->mParts[i]->id;
        nlohmann::json part = p->mParts[i];
        partsDesc.push_back(part);
    }
    j["partOrder"] = partOrder;

    j["mParts"] = partsDesc;

    //Then components and bonds
    j["mComponents"] = p->mComponents;
    std::map<int, std::string> componentOrder;
    for (int i = 0; i < p->mComponents.size(); i++) {
        componentOrder[i] = p->mComponents[i]->getId();
    }
    j["componentOrder"] = componentOrder;
    j["mBonds"] = p->mBonds;
    //Make a note if state equations have been computed (sometimes unintialized densmatrix will have nonzero rows but no columns)
    j["computed"] = p->computedCoordinates.nrows() > 0 && p->computedCoordinates.ncols() == 1;
}

void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraph> &p)
{
    //First do the child bond graphs
    if (j["class"] == "BondGraph") {
        j["id"].get_to(p->id);
        j["cloned"].get_to(p->cloned);
        std::vector<RCPLIB::RCP<BondGraph>> mParts;

        if (j["mParts"].size() > 0) {
            //Load subparts
            std::map<std::string, RCPLIB::RCP<BondGraph>> parts;
            //Ensure we are loading in the correct order
            //Json does not maintian ordered sets
            std::map < std::string, nlohmann::json> partJSONs;
            for (int i = 0; i < j["mParts"].size(); i++)
            {
                partJSONs[(j["mParts"][i])["id"]] = (j["mParts"][i]);
            }
            std::map<int, std::string> partOrder;
            j["partOrder"].get_to(partOrder);
            for (int i = 0; i < partOrder.size(); i++) {
                RCPLIB::RCP<BondGraph> part = createBondGraph();
                //Remove from registry as id will change, addBondGraph will register with new id                
                ComponentRegistry::defaultRegistry().removeComponent(part->id);
                partJSONs[partOrder[i]].get_to(part);
                p->addBondGraph(part);
                //Set owner post addition else addBondGraph will complain
                j["ownerID"].get_to(p->ownerID);
            }
        }
        //Add the parts
        std::unordered_map<std::string, RCPLIB::RCP<BondGraphElementBase>> comps;
        std::unordered_map<std::string, RCPLIB::RCP<BondGraphElementBase>> proxy;
        std::unordered_map<std::string, bool> flags;
        int cnt = 0;
        for (auto &c : p->mComponents) {
            cnt++;
            if (c->isProxy()) {
                proxy[c->getId()] = c;
                flags[c->getId()] = true;
            } else {
                comps[c->getId()] = c;
                flags[c->getId()] = false;
            }
        }
        

        //Then components and bonds
        bool proxyFlag;
        for (auto c : j["mComponents"]) {
            c["proxy"].get_to(proxyFlag);
            if (!proxyFlag && comps.find(c["data"]["mId"])==comps.end()) { //Create if it is not already there
                const RCPLIB::RCP<BondGraphElementBase> &comp = loadJson(c, RCPLIB::null);
                comps[comp->getId()] = comp;
                flags[comp->getId()] = false;
            }
        }
        for (auto c : j["mComponents"]) {
            c["proxy"].get_to(proxyFlag);
            if (proxyFlag && proxy.find(c["data"]["mId"])==proxy.end()) {//Create if it is not already in there
                const RCPLIB::RCP<BondGraphElementBase> &pro = comps[c["data"]["mId"]];
                assert(pro.is_null()==false);
                const RCPLIB::RCP<BondGraphElementBase> &comp = loadJson(c, pro->_data);
                proxy[comp->getId()] = comp;
                flags[comp->getId()] = true;
            }
        }

        //Get the component ids to create bonds
        bool pflag;
        std::map<std::string, RCPLIB::RCP<BondGraphElementBase>> cMap;
        std::map<std::string,bool> existingBonds;
        for(auto b : p->mBonds){
            existingBonds[b->bid] = true;
        }
        for (auto b : j["mBonds"]) {
            /**
             * @brief Proxy handling
             * Both the core and the proxy have the same port vector, so its immaterial which
             * handle is obtained as long as it exists
             */
            auto fromComp = comps[b["mFromPort"]["mComponent"]];
            if (fromComp.is_null() && ! proxy[b["mFromPort"]["mComponent"]].is_null()){
                fromComp = proxy[b["mFromPort"]["mComponent"]];
            }
            int fpin = b["mFromPort"]["portNo"];
            auto toComp = comps[b["mToPort"]["mComponent"]];
            if (toComp.is_null() && !proxy[b["mToPort"]["mComponent"]].is_null())
                toComp = proxy[b["mToPort"]["mComponent"]];
                
            int tpin = b["mToPort"]["portNo"];

            if (cMap.find(fromComp->mId) == cMap.end() && fromComp->isProxy() == flags[fromComp->mId])
                cMap[fromComp->getId()] = fromComp;
            if (cMap.find(toComp->mId) == cMap.end() && toComp->isProxy() == flags[toComp->mId])
                cMap[toComp->getId()] = toComp;

            //Create if not already present
            if(existingBonds[b["bid"]]){
                continue;
            }
            auto newBond = createBond(b["mId"]);
            b["mName"].get_to(newBond->mName);
            b["mStorageType"].get_to(newBond->mStorageType);
            b["mGroup"].get_to(newBond->mGroup);
            b["mStorageIndex"].get_to(newBond->mStorageIndex);
            b["bid"].get_to(newBond->bid);

            newBond->connect(fromComp->mPorts[fpin], toComp->mPorts[tpin]);
            p->addBond(newBond);                
        }
        //Load the used component handles
        std::map<int, std::string> componentOrder;
        j["componentOrder"].get_to(componentOrder);
        for (int i = 0; i < componentOrder.size();i++) {
            //Only add of not already found
            auto co = cMap[componentOrder[i]];
            if(std::find(p->mComponents.begin(),p->mComponents.end(),co)==p->mComponents.end()){
                p->mComponents.push_back(co);
            }
        }
        //Check if the state equations were computed, if so compute
        j["computed"].get_to(pflag);
        if (pflag)
            p->computeStateEquation();
        std::string readableName = j["readableName"];
        if (readableName != "")
            ComponentRegistry::defaultRegistry().setName(p, readableName);
        ComponentRegistry::defaultRegistry().ownBondgraph(j["ownerID"], p, readableName);
    } else {
        throw BGException("Invalid JSON");
    }
}

void from_json(const nlohmann::json &j, units::precise_unit &p)
{
    p = units::unit_from_string(j);
}

void to_json(nlohmann::json &j, const units::precise_unit &p)
{
    j = units::to_string(p);
}

void from_json(const nlohmann::json &j, RCPLIB::RCP<Value> &p)
{
    j["prefix"].get_to(p->prefix);
    j["name"].get_to(p->name);
    j["value"].get_to(p->value);
    //j["units"].get_to(p->units);
    if (j["units"] != "")
        p->units = units::unit_from_string(j["units"]);
    else
        p->units = units::precise::one;
    j["universalConstant"].get_to(p->universalConstant);
}

RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &proxy)
{
    if (j["class"] != "BondGraphElementBase")
        throw BGException("Invalid JSON");
    RCPLIB::RCP<ElementImpl> data = proxy;
    if (j["proxy"] && proxy.is_null()) {
        throw BGException("JSON belongs to a proxy, but pointer to candidate data not provided");
    }

    RCPLIB::RCP<BondGraphElementBase> ptr;
    const RCPLIB::RCP<ElementImpl>& obj = proxy; //Lambda function requires obj
    PassiveType etype;
    j["data"]["mElementType"].get_to(etype);
    //Define
    #include "createElements.func"
    if(etype==ePortHamiltonian){
        const RCPLIB::RCP<PortHamiltonian> &ph = RCPLIB::rcp_dynamic_cast<PortHamiltonian>(ptr);
        j["constraints"].get_to(ph->constraints);
    }

    if (!j["proxy"]) {
        //Component with generated id is registered in the component registry
        //Deregister and register with given ide
        ComponentRegistry::defaultRegistry().removeComponent(ptr->mId);
        auto dj = j["data"];
        dj["mName"].get_to(ptr->mName);
        dj["mId"].get_to(ptr->mId);
        dj["hamiltonian"].get_to(ptr->hamiltonian);
        dj["mElementType"].get_to(ptr->mElementType);
        dj["mComponentGroup"].get_to(ptr->mComponentGroup);
        dj["mDomain"].get_to(ptr->mDomain);
        dj["numStates"].get_to(ptr->numStates);
        dj["allocatedPorts"].get_to(ptr->allocatedPorts);
        dj["portModified"].get_to(ptr->portModified);
        //Ports will be created based on type but set the names and properties based on these values
        dj["constitutiveEq"].get_to(ptr->constitutiveEq);
        dj["constitutiveEqIndex"].get_to(ptr->constitutiveEqIndex);
        auto params = dj["mParameter"];
        for (int i = 0; i < params.size(); i++) {
            std::string ix = std::to_string(i);
            params[ix].get_to(ptr->mParameter[i]);
        }

        const nlohmann::json &ports = dj["mPorts"];
        for (int i = 0; i < ports.size(); i++) {
            std::string ix = std::to_string(i);
            RCPLIB::RCP<Port> pv;
            if (j["data"]["allocatedPorts"]) {
                pv = ptr->mPorts[i];
                pv->connect(ptr.create_weak());
            } else {
                pv = createPort(ports[ix]["receivesPower"]);
                ptr->connect(pv); // Connect the component to the port
            }
            ports[ix]["receivesPower"].get_to(pv->receivesPower);
            ports[ix]["mdofIndex"].get_to(pv->mdofIndex);
            ports[ix]["weight"].get_to(pv->weight);
            ports[ix]["portNo"].get_to(pv->portNo);
        }
        std::string readableName = j["readableName"];
        ComponentRegistry::defaultRegistry().ownComponent(j["parent"], ptr, readableName);
    }
    j["parent"].get_to(ptr->parent);
    j["proxy"].get_to(ptr->proxy);
    j["dofID"].get_to(ptr->dofID);
    if (j["readableName"] != "")
        ComponentRegistry::defaultRegistry().setName(ptr->mId, j["readableName"]);
    logDebug("Loaded element ", ptr->mName, " type ", ptr->mElementType);
    return ptr;
}

} // namespace BG