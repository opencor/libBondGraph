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

/**
 * @brief Webinterface definition
 *
 * Functions to load bondgraphs descriptions in stringified json format and
 * generate cellml files are defined.
 */
#include <emscripten.h>
#include <emscripten/bind.h>
using namespace emscripten;
#include "bondgraph.hpp"
#include "thirdparty/json.hpp"

using namespace BG;

//Helper function to log progress in console
EM_JS(void, printInput, (const char *str), { console.log(UTF8ToString(str)); });


static RCPLIB::RCP<BondGraphInterface> generateBondGraph(std::string bgJson){
    nlohmann::json jf = nlohmann::json::parse(bgJson);
    nlohmann::json methods = getSupportedPhysicalDomainsAndFactoryMethods();
    std::map<std::string, RCPLIB::RCP<BG::BGElement>> elements;
    std::ostringstream ss;
    std::map<std::string, std::string> import_file;
    nlohmann::json items = jf["sceneItems"];
    std::map<std::string, std::string> importName;
    std::map<std::string, std::map<std::string, std::string>> sharedImports;
    std::map<std::string,
             std::map<std::string, RCPLIB::RCP<BG::BondGraphInterface>>>
        bondgraphInstances;
    std::map<std::string, RCPLIB::RCP<BG::BGElement>> proxies;
    std::vector<std::string> missingElements;

    auto ioBondGraph = createBondGraph();
    for (const auto &itm : items.items()) {
      auto k = itm.key();
      auto val = itm.value();
      std::string assignedImportName = "";
      if (importName.find(k) != importName.end()) {
        assignedImportName = importName[k];
      }
      if (val["typeId"] == "BGElement") {
        auto def = val["definition"];
        auto edef = def["annotation"]["ElementDefinition"];
        auto sp = def["annotation"]["statesAndParameters"];
        auto annot = def["annotation"]["Annotation"];
        std::string displayName = def["displayname"];
        //std::replace(displayName.begin(), displayName.end(), ':', 'c');

        auto cellmlName = displayName;
        auto dName = def["mid"]; // Connections store mid's

        std::string domain = edef["domain"];
        if(domain == "Annotation"){
          continue;
        }        
        std::string type = edef["type"];
        std::string clas = edef["class"];
        std::string mName = methods[domain][type];

        RCPLIB::RCP<BG::BGElement> bge;
        if (proxies.find(dName) != proxies.end()) {
          bge = proxies[dName];
        }
        if (clas == "userdefined") {
          auto mapDef = def["annotation"]["mapdefinition"];

          nlohmann::json vmap;
          vmap["type"] = "file";
          vmap["compartment"] = mapDef["component"];
          // Cut to basename
          std::string fileName = mapDef["uri"];
          std::string baseFilename =
              fileName.substr(fileName.find_last_of("/\\") + 1);
          import_file[baseFilename] = fileName;
          vmap["filename"] = baseFilename;
          vmap["target"] = mapDef["variable"];
          vmap["link"] = mapDef["link"];
          // If link is true provide "importName" name so that the same cellml
          // instance is linked
          //||assignedImportName!="" part is to ensure the main varible uses
          // the
          // same importName
          if (vmap["link"] || assignedImportName != "") {
            vmap["importName"] = assignedImportName;
          }
          if (mapDef.contains("timeVariable")) {
            vmap["mapvariables"] = mapDef["timeVariable"];
          }
          // Check if any varibles need to be mapped
          if (mapDef["sourceType"]) { // Potential type
            if (bge.is_null())
              bge = BG::createPotentialSource();
            bge->setParameter("u", vmap.dump(), mapDef["dimension"]);
            bge->setName(cellmlName);
          } else {
            if (bge.is_null())
              bge = BG::createFlowSource();
            bge->setName(cellmlName);
            bge->setParameter("i", vmap.dump(), mapDef["dimension"]);
          }
        } else if (clas == "junction") { // Handle TF and GY
          if (bge.is_null())
            bge = BG::createBondgraphElement(mName);
          bge->setName(cellmlName);
          if (sp.is_object() || sp.is_array()) {
            for (const auto &pv : sp) {
              auto nm = pv["name"];
              if (pv["value"].is_object()) {
                auto pm = pv["value"];
                if (pm.contains("file")) { // When a file reference is made
                  nlohmann::json vmap;
                  vmap["type"] = "file";
                  vmap["compartment"] = pm["component"];
                  std::string fileName = pm["filename"];
                  std::string baseFilename =
                      fileName.substr(fileName.find_last_of("/\\") + 1);
                  vmap["filename"] = baseFilename;
                  vmap["target"] = pm["variable"];
                  vmap["link"] = pm["link"];
                  // If link is true provide "importName" name so that the
                  // same cellml instance is linked
                  if (vmap["link"] || assignedImportName != "") {
                    vmap["importName"] = assignedImportName;
                  }
                  if (pm.contains("timeVariable")) {
                    vmap["mapvariables"] = pm["timeVariable"];
                  }
                  // file, use for resolving multiple/shared instances
                  bge->setParameter(nm, vmap.dump(), "");
                } else { // User has changed just the value
                  ss.str("");
                  ss << pv["value"]["value"];
                  bge->setParameter(nm, ss.str(), "");
                }
              } else {
                ss.str("");
                ss << pv["value"];
                bge->setParameter(nm, ss.str(), "");
              }
            }
          }
        } else if (clas == "passive") { // If passive set prescribed state and
          // parameter values
          if (bge.is_null())
            bge = BG::createBondgraphElement(mName);
          bge->setName(cellmlName);

          for (const auto &pv : sp) {
            auto nm = pv["name"];
            if (pv["value"].is_object()) {
              auto pm = pv["value"];
              if (pm.contains("file")) { // When a file reference is made
                nlohmann::json vmap;
                vmap["type"] = "file";
                vmap["compartment"] = pm["component"];
                std::string fileName = pm["file"];
                std::string baseFilename =
                    fileName.substr(fileName.find_last_of("/\\") + 1);
                vmap["filename"] = baseFilename;
                vmap["target"] = pm["variable"];
                // Check if any varibles need to be mapped
                if (pm.contains("mapvariables")) {
                  vmap["mapvariables"] = pm["mapvariables"]["timeVariable"];
                }
                vmap["link"] = pm["link"];
                // If link is true provide "importName" name so that the same
                // cellml instance is linked
                if (vmap["link"] || assignedImportName != "") {
                  vmap["importName"] = assignedImportName;
                }
                if (pm.contains("timeVariable")) {
                  vmap["mapvariables"] = pm["timeVariable"];
                }
                // file, use for resolving multiple/shared instances
                bge->setParameter(nm, vmap.dump(), pv["dimension"]);
              } else { // User has changed just the value
                ss.str("");
                ss << pv["value"]["value"];
                bge->setParameter(nm, ss.str(), pv["dimension"]);
              }
            } else {
              ss.str("");
              ss << pv["value"];
              bge->setParameter(nm, ss.str(), pv["dimension"]);
            }
          }
        }
        bge->setPMRAnnotation(annot);
        elements[dName] = bge;
        ioBondGraph->addComponent(bge);
      }
    }
    // Create the connections
    for (const auto &itm : items.items()) {
      auto k = itm.key();
      auto val = itm.value();
      if (val["typeId"] == "BGDConnection") {
        auto def = val["definition"];
        auto fElem = elements[def["first"]];
        auto fPort = def["firstPort"];
        auto fports = fElem->getPorts();
        int cFP = 0;
        for (size_t pc = 0; pc < fports.size(); pc++) {
          if (fports[pc]->getId() == fPort) {
            cFP = (int)pc;
            break;
          }
        }
        auto lElem = elements[def["second"]];
        auto lPort = def["secondPort"];
        auto lports = lElem->getPorts();
        int cLP = 0;
        for (size_t pc = 0; pc < lports.size(); pc++) {
          if (lports[pc]->getId() == lPort) {
            cLP = (int)pc;
            break;
          }
        }
        auto fType = fElem->getType();
        auto lType = lElem->getType();
        bool fElemISReaction = fType == BG::PassiveType::bReaction;
        bool lElemISReaction = lType == BG::PassiveType::bReaction;
        bool fElemISJunction =
            fType == BG::PassiveType::eOne || fType == BG::PassiveType::eZero;
        bool lElemISJunction =
            lType == BG::PassiveType::eOne || lType == BG::PassiveType::eZero;

        if (!fElemISReaction && !lElemISReaction) {
          if (!fElemISJunction && !lElemISJunction) {
            ioBondGraph->connect(fElem, cFP, lElem, cLP);
          } else if (fElemISJunction && !lElemISJunction) {
            ioBondGraph->connect(
                fElem, lElem, 0,
                cLP); // Match the signature to avoid the ambiguity
          } else if (!fElemISJunction && lElemISJunction) {
            ioBondGraph->connect(fElem, cFP, lElem);
          } else {
            ioBondGraph->connect(fElem, lElem);
          }
        } else if (fElemISReaction) {
          ioBondGraph->connectInverting(fElem, 1, lElem);
        } else {
          ioBondGraph->connectInverting(fElem, lElem, 0);
        }
      }
    }
  return ioBondGraph;
}

std::string getSupportedPhysicalDomainsAndFactories(){
  nlohmann::json methods = getSupportedPhysicalDomainsAndFactoryMethods();
  std::string result = methods.dump();
  return result;
}

std::string generateCellMLForBondgraph(std::string bgJson) {
  nlohmann::json result;
  result["success"] = true;
  result["input"] = bgJson;

  try{
    auto ioBondGraph = generateBondGraph(bgJson);
    printInput("Successfully generated BondGraph!");
    auto eqs = ioBondGraph->computeStateEquation();
    printInput("Successfully computed state equation!");
    //Get project name
    nlohmann::json jf = nlohmann::json::parse(bgJson);
    if(jf["Provenance"].contains("projectname")){
      std::string projectName = jf["Provenance"]["projectname"];
      auto files = getCellML(projectName, ioBondGraph, eqs);
      result["cellml"] = files;
      printInput("Successfully generated cellml!");
    }else{
      auto files = getCellML("WASM", ioBondGraph, eqs);
      result["cellml"] = files;
      printInput("Successfully generated cellml with temporary projectname `WASM`!");
    }
    
    
  } catch (const std::exception &exc) {
    result["success"] = false;
    result["error"] = exc.what();
    printInput(exc.what());
  } catch (nlohmann::json::parse_error &ex) {
    result["success"] = false;
    result["error"] = ex.what();
    printInput(ex.what());
  }

  std::string output = result.dump();
  return output;
}

std::string generatePortHamiltonian(std::string bgJson){
  nlohmann::json result;
  result["success"] = true;
  result["input"] = bgJson;

  try{
    auto ioBondGraph = generateBondGraph(bgJson);
    printInput("Successfully generated BondGraph!");
    auto phs = ioBondGraph->computePortHamiltonian();
    printInput("Successfully computed portHamiltonian!");
    result["phs"] = phs;
  } catch (const std::exception &exc) {
    result["success"] = false;
    result["error"] = exc.what();
    printInput(exc.what());
  } catch (nlohmann::json::parse_error &ex) {
    result["success"] = false;
    result["error"] = ex.what();
    printInput(ex.what());
  }
  std::string output = result.dump();
  return output;
}

EMSCRIPTEN_BINDINGS(libbondgraph) {
  function("getSupportedPhysicalDomainsAndFactoryMethods",&getSupportedPhysicalDomainsAndFactories);
  function("generateCellMLForBondgraph", &generateCellMLForBondgraph);
  function("generatePortHamiltonian",&generatePortHamiltonian);
  function("checkUnits",&BG::checkUnits);
  //To create Bondgraphs - send the scene and create it in the library
  //and send the json
}