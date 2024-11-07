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

// Helper function to log progress in console
EM_JS(void, printInput, (const char *str), { console.log(UTF8ToString(str)); });

std::string getSupportedPhysicalDomainsAndFactories() {
  nlohmann::json methods = BG::getSupportedPhysicalDomainsAndFactoryMethods();
  std::string result = methods.dump();
  return result;
}

std::string generateCellMLForBondgraph(std::string bgJson) {
  nlohmann::json result;
  result["success"] = true;
  result["input"] = bgJson;

  try {
    auto ioBondGraph = BG::generateBondGraph(bgJson);
    printInput("Successfully generated BondGraph!");
    auto eqs = ioBondGraph->computeStateEquation();
    printInput("Successfully computed state equation!");
    // Get project name
    nlohmann::json jf = nlohmann::json::parse(bgJson);
    if (jf["Provenance"].contains("projectname")) {
      std::string projectName = jf["Provenance"]["projectname"];
      auto files = BG::getCellML(projectName, ioBondGraph, eqs);
      result["cellml"] = files;
      printInput("Successfully generated cellml!");
    } else {
      auto files = BG::getCellML("WASM", ioBondGraph, eqs);
      result["cellml"] = files;
      printInput(
          "Successfully generated cellml with temporary projectname `WASM`!");
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

std::string generatePortHamiltonian(std::string bgJson) {
  nlohmann::json result = BG::generatePortHamiltonian(bgJson);
  if (result.contains("error")) {
    std::string err = result["error"];
    printInput(err.c_str());
  }
  std::string output = result.dump();
  return output;
}

EMSCRIPTEN_BINDINGS(libbondgraph) {
  function("getSupportedPhysicalDomainsAndFactoryMethods",
           &getSupportedPhysicalDomainsAndFactories);
  function("generateCellMLForBondgraph", &generateCellMLForBondgraph);
  function("generatePortHamiltonian", &generatePortHamiltonian);
  function("checkUnits", &BG::checkUnits);
  // To create Bondgraphs - send the scene and create it in the library
  // and send the json
}