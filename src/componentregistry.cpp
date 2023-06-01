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


#include "Elementsbase.h"
#include "Exceptions.h"
#include "bondgraph.h"
#include "cellmlunitsmap.h"
#include "componentregistry.h"
#include <stdlib.h> //For srand
#include <time.h>


namespace BG {

ComponentRegistryInterface &defaultRegistry()
{
    return ComponentRegistry::defaultRegistry();
}

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

void newWorkSpace()
{
    ComponentRegistry::defaultRegistry().newWorkSpace();
}

void ComponentRegistry::newWorkSpace()
{
    if (!instance.is_null()) {
        instance->clearWorkspace();
    }
}

void createRegistry()
{
    ComponentRegistry::defaultRegistry().createRegistry();
}

void ComponentRegistry::createRegistry()
{
    if (instance.is_null()) {
        instance = RCPLIB::rcp(new ComponentRegistry());
        logInfo("Registry created");
        instance->registryCleared = true;
    }
}

ComponentRegistry::ComponentRegistry()
{
    srand(time(NULL));
    clearWorkspace();
    //logCritical("Random number seed set to 0; modify this at line 49, componentregsitry.cpp");
    //srand(0);
}

ComponentRegistry::~ComponentRegistry()
{
    clearWorkspace();
}

void ComponentRegistry::ownComponent(std::string parent, RCPLIB::RCP<BGElement> comp, std::string readableName)
{
    comp.set_has_ownership();
    addComponent(parent, comp, readableName);
}

void ComponentRegistry::addComponent(std::string parent, const RCPLIB::RCP<BGElement> &comp, std::string readableName)
{
    if (instance->components.find(comp->getId()) == instance->components.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->components[comp->getId()] = std::make_tuple(parent, comp, json);
        if (readableName != "") {
            instance->nameMap[readableName] = comp->getId();
            instance->nameMap[comp->getId()] = readableName;
        }
    } else {
        for (auto x : instance->components) {
            auto elem = std::get<1>(x.second);
            logDebug("Existing Element ID ", x.first, " element ", elem->getId(), " displayname ", instance->nameMap[elem->getId()]);
        }
        throw BGException("Component with id " + comp->getId() + " already exists");
    }
}

void ComponentRegistry::ownBondgraph(std::string parent, RCPLIB::RCP<BondGraphInterface> comp, std::string readableName)
{
    comp.set_has_ownership();
    addBondgraph(parent, comp, readableName);
}

void ComponentRegistry::addBondgraph(std::string parent, const RCPLIB::RCP<BondGraphInterface> &comp, std::string readableName)
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
void ComponentRegistry::addOrReplaceComponent(std::string parent, RCPLIB::RCP<BGElement> comp, std::string readableName_)
{
    comp.set_has_ownership();
    std::string readableName = readableName_;
    if (instance->components.find(comp->getId()) == instance->components.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->components[comp->getId()] = std::make_tuple(parent, comp, json);
    } else {
        std::string par;
        RCPLIB::RCP<BGElement> com;
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
void ComponentRegistry::addOrReplaceBondgraph(std::string parent, RCPLIB::RCP<BondGraphInterface> comp, std::string readableName_)
{
    comp.set_has_ownership();
    std::string readableName = readableName_;
    if (instance->bondgraphs.find(comp->getId()) == instance->bondgraphs.end()) {
        RCPLIB::RCP<nlohmann::json> json = RCPLIB::rcp(new nlohmann::json);
        instance->bondgraphs[comp->getId()] = std::make_tuple(parent, comp, json);
    } else {
        std::string par;
        RCPLIB::RCP<BondGraphInterface> com;
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
void ComponentRegistry::setName(const RCPLIB::RCP<BGElement> &comp, std::string readableName)
{
    setName(comp->getId(), readableName);
}
void ComponentRegistry::setName(const RCPLIB::RCP<BondGraphInterface> &comp, std::string readableName)
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
std::tuple<int, const RCPLIB::RCP<BGElement> &, const RCPLIB::RCP<BondGraphInterface> &> ComponentRegistry::getComponentByName(std::string readableName)
{
    if (instance->nameMap.find(readableName) != instance->nameMap.end()) {
        std::string id = nameMap[readableName];
        if (instance->components.find(id) != instance->components.end()) {
            RCPLIB::RCP<BondGraphInterface> ptr;
            return std::make_tuple(0, std::get<1>(instance->components[id]), ptr);
        } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
            RCPLIB::RCP<BGElement> ptr;
            return std::make_tuple(1, ptr, std::get<1>(instance->bondgraphs[id]));
        } else {
            throw BGException("Component with readableName " + readableName + " and id " + id + " doesn't exist!");
        }
    } else {
        throw BGException("Entry with readableName " + readableName + " doesn't exist!");
    }
}
std::tuple<int, const RCPLIB::RCP<BGElement> &, const RCPLIB::RCP<BondGraphInterface> &> ComponentRegistry::getComponent(std::string id)
{
    if (instance->components.find(id) != instance->components.end()) {
        RCPLIB::RCP<BondGraphInterface> ptr;
        return std::make_tuple(0, std::get<1>(instance->components[id]), ptr);
    } else if (instance->bondgraphs.find(id) != instance->bondgraphs.end()) {
        RCPLIB::RCP<BGElement> ptr;
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

std::unordered_map<std::string, std::string> &ComponentRegistry::getNameMap()
{
    return nameMap;
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

RegistryCreator createRegistryInstance;

} // namespace BG