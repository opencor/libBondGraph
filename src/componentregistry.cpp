#include "Exceptions.h"
#include "cellmlunitsmap.h"
#include "componentregistry.h"
#include <stdlib.h> //For srand
#include <time.h>

namespace BG {
ComponentRegistry &ComponentRegistry::defaultRegistry()
{
    return *instance;
}

void ComponentRegistry::clearWorkspace()
{
    for (auto &c : components) {
        std::get<1>(c.second) = RCPLIB::null;
        std::get<2>(c.second) = RCPLIB::null;
    }
    for (auto &c : bondgraphs) {
        std::get<1>(c.second) = RCPLIB::null;
        std::get<2>(c.second) = RCPLIB::null;
    }

    nameMap.clear();
    components.clear();
    bondgraphs.clear();
}

void ComponentRegistry::newWorkSpace()
{
    if (!instance.is_null()) {
        instance->clearWorkspace();
        auto curMem = mallinfo();
        logInfo("Registry cleared, freeing ", curMem.uordblks - instance->memUsage.uordblks, " bytes of memory");
        instance->memUsage = mallinfo();
    }
}

void ComponentRegistry::createRegistry()
{
    if (instance.is_null()) {
        instance = RCPLIB::rcp(new ComponentRegistry());
        logInfo("Registry created");
        instance->memUsage = mallinfo();
        instance->registryCleared = true;
    }
}

ComponentRegistry::ComponentRegistry()
{
    srand(time(NULL));
    //logCritical("Random number seed set to 0; modify this at line 49, componentregsitry.cpp");
    //srand(0);
}

ComponentRegistry::~ComponentRegistry()
{
    clearWorkspace();
}

void ComponentRegistry::ownComponent(std::string parent, RCPLIB::RCP<BondGraphElementBase> comp, std::string readableName)
{
    comp.set_has_ownership();
    addComponent(parent, comp, readableName);
}

void ComponentRegistry::addComponent(std::string parent, const RCPLIB::RCP<BondGraphElementBase> &comp, std::string readableName)
{
    if (instance->components.find(comp->getId()) == instance->components.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->components[comp->getId()] = std::make_tuple(parent, comp, json);
        if (readableName != "") {
            instance->nameMap[readableName] = comp->getId();
            instance->nameMap[comp->getId()] = readableName;
        }
    } else {
        throw BGException("Component with id " + comp->getId() + " already exists");
    }
}

void ComponentRegistry::ownBondgraph(std::string parent, RCPLIB::RCP<BondGraph> comp, std::string readableName)
{
    comp.set_has_ownership();
    addBondgraph(parent, comp, readableName);
}

void ComponentRegistry::addBondgraph(std::string parent, const RCPLIB::RCP<BondGraph> &comp, std::string readableName)
{
    if (instance->bondgraphs.find(comp->getId()) == instance->bondgraphs.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->bondgraphs[comp->getId()] = std::make_tuple(parent, comp, json);
        if (readableName != "") {
            instance->nameMap[readableName] = comp->getId();
            instance->nameMap[comp->getId()] = readableName;
        }
    } else {
        throw BGException("Bondgraph with id " + comp->getId() + " already exists");
    }
}
void ComponentRegistry::addOrReplaceComponent(std::string parent, RCPLIB::RCP<BondGraphElementBase> comp, std::string readableName_)
{
    comp.set_has_ownership();
    std::string readableName = readableName_;
    if (instance->components.find(comp->getId()) == instance->components.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->components[comp->getId()] = std::make_tuple(parent, comp, json);
    } else {
        std::string par;
        RCPLIB::RCP<BondGraphElementBase> com;
        RCPLIB::RCP<nlohmann::json> json;
        std::tie(par, com, json) = instance->components[comp->getId()];
        instance->components.erase(comp->getId());
        instance->components[comp->getId()] = std::make_tuple(parent, comp, json);
        if (readableName_ == "") {
            if (instance->nameMap.find(comp->getId()) != instance->nameMap.end()) {
                readableName = instance->nameMap[comp->getId()];
            }
        }
    }
    if (readableName != "") {
        instance->nameMap[readableName] = comp->getId();
        instance->nameMap[comp->getId()] = readableName;
    }
}
void ComponentRegistry::addOrReplaceBondgraph(std::string parent, RCPLIB::RCP<BondGraph> comp, std::string readableName_)
{
    comp.set_has_ownership();
    std::string readableName = readableName_;
    if (instance->bondgraphs.find(comp->getId()) == instance->bondgraphs.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->bondgraphs[comp->getId()] = std::make_tuple(parent, comp, json);
    } else {
        std::string par;
        RCPLIB::RCP<BondGraph> com;
        RCPLIB::RCP<nlohmann::json> json;
        std::tie(par, com, json) = instance->bondgraphs[comp->getId()];
        instance->bondgraphs.erase(comp->getId());
        instance->bondgraphs[comp->getId()] = std::make_tuple(parent, comp, json);
        if (readableName_ == "") {
            if (instance->nameMap.find(comp->getId()) != instance->nameMap.end()) {
                readableName = instance->nameMap[comp->getId()];
            }
        }
    }
    if (readableName != "") {
        instance->nameMap[readableName] = comp->getId();
        instance->nameMap[comp->getId()] = readableName;
    }
}

void ComponentRegistry::removeComponentByName(std::string readableName)
{
    if (instance->nameMap.find(readableName) != instance->nameMap.end()) {
        std::string id = nameMap[readableName];
        if (instance->components.find(id) != instance->components.end()) {
            instance->components.erase(id);
            instance->nameMap.erase(readableName);
            instance->nameMap.erase(id);
        } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
            instance->bondgraphs.erase(id);
            instance->nameMap.erase(readableName);
            instance->nameMap.erase(id);
        } else {
            instance->nameMap.erase(readableName);
            throw BGException("Component with readableName " + readableName + " and id " + id + " doesn't exist!");
        }
    } else {
        throw BGException("Entry with readableName " + readableName + " doesn't exist!");
    }
}

void ComponentRegistry::removeComponentIgnoreIfNotAvailable(std::string id)
{
    if (instance->components.find(id) != instance->components.end()) {
        instance->components.erase(id);
        auto readableName = nameMap[id];
        instance->nameMap.erase(readableName);
        instance->nameMap.erase(id);
    } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
        instance->bondgraphs.erase(id);
        auto readableName = nameMap[id];
        instance->nameMap.erase(readableName);
        instance->nameMap.erase(id);
    }
}

void ComponentRegistry::removeComponent(std::string id)
{
    if (instance->components.find(id) != instance->components.end()) {
        instance->components.erase(id);
        auto readableName = nameMap[id];
        instance->nameMap.erase(readableName);
        instance->nameMap.erase(id);
    } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
        instance->bondgraphs.erase(id);
        auto readableName = nameMap[id];
        instance->nameMap.erase(readableName);
        instance->nameMap.erase(id);
    } else {
        auto readableName = nameMap[id];
        instance->nameMap.erase(readableName);
        instance->nameMap.erase(id);
        throw BGException("Component with readableName " + readableName + " and id " + id + " doesn't exist!");
    }
}
void ComponentRegistry::setName(const RCPLIB::RCP<BondGraphElementBase> &comp, std::string readableName)
{
    setName(comp->getId(), readableName);
}
void ComponentRegistry::setName(const RCPLIB::RCP<BondGraph> &comp, std::string readableName)
{
    setName(comp->getId(), readableName);
}

void ComponentRegistry::setName(std::string id, std::string readableName)
{
    if (instance->nameMap.find(id) != instance->nameMap.end()) {
        auto readableName = nameMap[id];
        instance->nameMap.erase(readableName);
        instance->nameMap.erase(id);
        instance->nameMap[id] = readableName;
        instance->nameMap[readableName] = id;
    }
}
std::tuple<int, const RCPLIB::RCP<BondGraphElementBase> &, const RCPLIB::RCP<BondGraph> &> ComponentRegistry::getComponentByName(std::string readableName)
{
    if (instance->nameMap.find(readableName) != instance->nameMap.end()) {
        std::string id = nameMap[readableName];
        if (instance->components.find(id) != instance->components.end()) {
            RCPLIB::RCP<BondGraph> ptr;
            return std::make_tuple(0, std::get<1>(instance->components[id]), ptr);
        } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
            RCPLIB::RCP<BondGraphElementBase> ptr;
            return std::make_tuple(1, ptr, std::get<1>(instance->bondgraphs[id]));
        } else {
            throw BGException("Component with readableName " + readableName + " and id " + id + " doesn't exist!");
        }
    } else {
        throw BGException("Entry with readableName " + readableName + " doesn't exist!");
    }
}
std::tuple<int, const RCPLIB::RCP<BondGraphElementBase> &, const RCPLIB::RCP<BondGraph> &> ComponentRegistry::getComponent(std::string id)
{
    if (instance->components.find(id) != instance->components.end()) {
        RCPLIB::RCP<BondGraph> ptr;
        return std::make_tuple(0, std::get<1>(instance->components[id]), ptr);
    } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
        RCPLIB::RCP<BondGraphElementBase> ptr;
        return std::make_tuple(1, ptr, std::get<1>(instance->bondgraphs[id]));
    } else {
        throw BGException("Component with id " + id + " doesn't exist!");
    }
}
RCPLIB::RCP<nlohmann::json> ComponentRegistry::getPropertiesByName(std::string readableName)
{
    if (instance->nameMap.find(readableName) != instance->nameMap.end()) {
        std::string id = nameMap[readableName];
        if (instance->components.find(id) != instance->components.end()) {
            return std::get<2>(instance->components[id]);
        } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
            return std::get<2>(instance->bondgraphs[id]);
        } else {
            throw BGException("Component with readableName " + readableName + " and id " + id + " doesn't exist!");
        }
    } else {
        throw BGException("Entry with readableName " + readableName + " doesn't exist!");
    }
}
RCPLIB::RCP<nlohmann::json> ComponentRegistry::getProperties(std::string id)
{
    if (instance->components.find(id) != instance->components.end()) {
        return std::get<2>(instance->components[id]);
    } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
        return std::get<2>(instance->bondgraphs[id]);
    } else {
        throw BGException("Component with id " + id + " doesn't exist!");
    }
}

RCPLIB::RCP<ComponentRegistry> ComponentRegistry::instance = RCPLIB::null;
//RCPLIB::RCP<CellMLUnits> CellMLUnits::instance = RCPLIB::null;

//create an instance
class RegistryCreator
{
public:
    RegistryCreator()
    {
        ComponentRegistry::createRegistry();
        //CellMLUnits::createUnitsMapper();
    }
};

RegistryCreator createRegistry;

} // namespace BG