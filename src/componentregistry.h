#ifndef COMPONENT_REGISTRY_H_
#define COMPONENT_REGISTRY_H_
#include "Elementsbase.h"
#include "RCP.h"
#include "bondgraph.h"
#include "thirdparty/json.hpp"
#include <malloc.h>
#include <tuple>
#include <unordered_map>
namespace BG {

class ComponentRegistry
{
private:
    static RCPLIB::RCP<ComponentRegistry> instance;
    nlohmann::json constitutiveEquations;
    std::unordered_map<std::string, bool> supportedPhysicalDomains;
    struct mallinfo memUsage;
    bool registryCleared = false;
protected:
    void clearWorkspace();
    //Map of object id -> host bondgraph id, handle to element, properties json
    std::unordered_map<std::string, std::tuple<std::string, RCPLIB::RCP<BondGraphElementBase>, RCPLIB::RCP<nlohmann::json>>>
        components;
    std::unordered_map<std::string, std::tuple<std::string, RCPLIB::RCP<BondGraph>, RCPLIB::RCP<nlohmann::json>>> bondgraphs;
    std::unordered_map<std::string, std::string> nameMap;
    void addOrReplaceComponent(std::string parent, RCPLIB::RCP<BondGraphElementBase> comp, std::string readableName = "");
    void addOrReplaceBondgraph(std::string parent, RCPLIB::RCP<BondGraph> comp, std::string readableName = "");
    ComponentRegistry();

public:
    virtual ~ComponentRegistry();
    //These functions take ownership of component memory management
    void ownComponent(std::string parent, RCPLIB::RCP<BondGraphElementBase> comp, std::string readableName = "");
    void ownBondgraph(std::string parent, RCPLIB::RCP<BondGraph> comp, std::string readableName = "");
    //These functions do not take ownership of component memory management
    void addComponent(std::string parent, const RCPLIB::RCP<BondGraphElementBase> &comp, std::string readableName = "");
    void addBondgraph(std::string parent, const RCPLIB::RCP<BondGraph> &comp, std::string readableName = "");
    void setName(const RCPLIB::RCP<BondGraphElementBase> &comp, std::string readableName);
    void setName(const RCPLIB::RCP<BondGraph> &comp, std::string readableName);
    void setName(std::string id, std::string readableName);
    void removeComponentByName(std::string readableName);
    void removeComponent(std::string id);
    void removeComponentIgnoreIfNotAvailable(std::string id);
    std::tuple<int, const RCPLIB::RCP<BondGraphElementBase> &, const RCPLIB::RCP<BondGraph> &> getComponentByName(std::string readableName);
    std::tuple<int, const RCPLIB::RCP<BondGraphElementBase> &, const RCPLIB::RCP<BondGraph> &> getComponent(std::string id);
    RCPLIB::RCP<nlohmann::json> getPropertiesByName(std::string readableName);
    RCPLIB::RCP<nlohmann::json> getProperties(std::string id);

    friend class BondGraph;
    static ComponentRegistry &defaultRegistry();
    static void createRegistry();
    /**
     * @brief Clears existing workspace, frees memory and starts afresh
     * User is responsible for saving any Changes in the objects prior to calling this function
     */
    static void newWorkSpace();

    std::unordered_map<std::string, std::string> &getNameMap()
    {
        return nameMap;
    }
};

} // namespace BG

#endif