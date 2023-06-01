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
#include <tuple>
#include <unordered_map>

namespace BG {

class ComponentRegistry : public ComponentRegistryInterface
{
private:
    static RCPLIB::RCP<ComponentRegistry> instance;
    nlohmann::json constitutiveEquations;
    std::unordered_map<std::string, bool> supportedPhysicalDomains;
    bool registryCleared = false;
protected:
    void clearWorkspace();
    //Map of object id -> host bondgraph id, handle to element, properties json
    std::unordered_map<std::string, std::tuple<std::string, RCPLIB::RCP<BGElement> , RCPLIB::RCP<nlohmann::json>>>
        components;
    std::unordered_map<std::string, std::tuple<std::string, RCPLIB::RCP<BondGraphInterface>, RCPLIB::RCP<nlohmann::json>>> bondgraphs;
    std::unordered_map<std::string, std::string> nameMap;
    void addOrReplaceComponent(std::string parent, RCPLIB::RCP<BGElement>  comp, std::string readableName = "");
    void addOrReplaceBondgraph(std::string parent, RCPLIB::RCP<BondGraphInterface> comp, std::string readableName = "");
    ComponentRegistry();

public:
    virtual ~ComponentRegistry();
    //These functions take ownership of component memory management
    void ownComponent(std::string parent, RCPLIB::RCP<BGElement>  comp, std::string readableName = "");
    void ownBondgraph(std::string parent, RCPLIB::RCP<BondGraphInterface> comp, std::string readableName = "");
    //These functions do not take ownership of component memory management
    void addComponent(std::string parent, const RCPLIB::RCP<BGElement>  &comp, std::string readableName = "");
    void addBondgraph(std::string parent, const RCPLIB::RCP<BondGraphInterface> &comp, std::string readableName = "");
    void setName(const RCPLIB::RCP<BGElement>  &comp, std::string readableName);
    void setName(const RCPLIB::RCP<BondGraphInterface> &comp, std::string readableName);
    void setName(std::string id, std::string readableName);
    void removeComponentByName(std::string readableName);
    void removeComponent(std::string id);
    void removeComponentIgnoreIfNotAvailable(std::string id);
    std::tuple<int, const RCPLIB::RCP<BGElement>  &, const RCPLIB::RCP<BondGraphInterface> &> getComponentByName(std::string readableName);
    std::tuple<int, const RCPLIB::RCP<BGElement>  &, const RCPLIB::RCP<BondGraphInterface> &> getComponent(std::string id);
    RCPLIB::RCP<nlohmann::json> getPropertiesByName(std::string readableName);
    RCPLIB::RCP<nlohmann::json> getProperties(std::string id);
    std::unordered_map<std::string, std::string> &getNameMap();
    
    friend class BondGraph;
    static ComponentRegistry &defaultRegistry();
    static void createRegistry();
    /**
     * @brief Clears existing workspace, frees memory and starts afresh
     * User is responsible for saving any Changes in the objects prior to calling this function
     */
    static void newWorkSpace();
};

} // namespace BG
