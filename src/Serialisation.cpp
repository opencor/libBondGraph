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

#include "Exceptions.h"
#include "bondgraph.h"
#include "cellmlunitsmap.h"
#include "thirdparty/tinyxml2.h"
#include <algorithm>
#include <sstream>
#include <symengine/eval.h>
#include <symengine/matrix.h>
#include <symengine/parser.h>
#include <symengine/parser/parser.h>
#include <symengine/printers.h>
#include <symengine/simplify.h>
#include <symengine/solve.h>
#include <symengine/visitor.h>
#include <tuple>

#include <boost/regex.hpp>

namespace BG {
std::string getMathML(const RCPLIB::RCP<const SymEngine::Basic> &expr) {
  std::string expression = SymEngine::mathml(*expr);
  const std::string from = "type=\"integer\"";
  const std::string to = "cellml:units=\"dimensionless\"";

  auto toCellMLUnit = [&expression, from, to]() {
    size_t start_pos = 0;
    while ((start_pos = expression.find(from, start_pos)) !=
           std::string::npos) {
      expression.replace(start_pos, from.length(), to);
      start_pos +=
          to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return expression;
  };

  toCellMLUnit();
  return expression;
}

std::vector<std::string> getMathML(SymEngine::vec_basic &expressions) {
  std::vector<std::string> mexpr;
  std::ostringstream ss;
  for (auto expr : expressions) {
    ss.str("");
    ss.clear();
    ss << "<apply><eq/>" << getMathML(expr)
       << "<cn cellml:units=\"dimensionless\">0.0</cn></apply>";
    std::string expression = ss.str();
    mexpr.push_back(expression);
  }
  return mexpr;
}

std::vector<std::string>
getMathML(std::map<RCPLIB::RCP<const SymEngine::Basic>,
                   RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
              &expressions) {
  std::vector<std::string> mexpr;
  std::ostringstream ss;
  for (auto expr : expressions) {
    ss.str("");
    ss.clear();
    ss << *expr.first;

    std::string term = ss.str();
    ss.str("");
    ss.clear();
    ss << "<apply><eq/>" << getMathML(expr.first) << getMathML(expr.second)
       << "</apply>";
    mexpr.push_back(ss.str());
  }
  return mexpr;
}

/**
 * @brief Generate the cellml code for import and linking external component
 *
 * @param vname Name of the variable that needs to be linked
 * @param prop Properties associated with that variable including the JSON
 * description of the target file import statement, connection statement, cellml
 * statement (if import statement, this will be ""), variable name, timevariable
 * name to be mapped (if none, this will be "") bool (true if dot_vname is to be
 * mapped)
 * @return std::tuple<std::string, std::string, std::string,
 * std::string,std::string,bool>
 */
std::tuple<std::string, std::string, std::string, std::string, std::string,
           bool>
generateLinkages(std::string vname,
                 std::tuple<std::string, std::string, char> &prop) {
  std::string value = std::get<1>(prop);
  // std::cerr << vname
  // <<"\t"<<std::get<0>(prop)<<"\t"<<std::get<1>(prop)<<"\t"<<std::get<2>(prop)<<std::endl;
  try {
    nlohmann::json j = nlohmann::json::parse(value);

    if (std::get<2>(prop) == 'c') { // Control variable
      // Has vname and dot_vname
      std::string compartment = j["compartment"];
      std::string targetVar = vname;
      if (j.contains("target")) {
        targetVar = j["target"];
      }
      std::string referenceName = vname;
      if (j.contains("importName")) {
        referenceName = j["importName"];
      }
      std::string filename = j["filename"];
      std::replace(filename.begin(), filename.end(), ' ', '_');
      std::string import = "<import xlink:href = \"" + filename + "\" >";
      import = import + "<component component_ref = \"" + compartment +
               "\" name = \"" + referenceName + "_IMP\" /> </import>";
      // std::string connection = "<connection>\n<map_components component_1 =
      // \"" + vname + "_IMP\" component_2 = \"main\" />"; connection =
      // connection + "<map_variables variable_1 = \""+targetVar+"\" variable_2
      // = \""+vname+"\" />";
      std::string connection = "<map_variables variable_1 = \"" + targetVar +
                               "\" variable_2 = \"" + vname + "\" />";
      if (j.contains("stateValue") &&
          j["stateValue"]) { // State variable mapping
        // in this case map to a variable named vname_smap, which will be
        // provided as init for the state variable
        connection = "<map_variables variable_1 = \"" + targetVar +
                     "\" variable_2 = \"" + vname + "_smap\" />";
      }

      bool derivativeMapped = false;
      std::string timeVariable = "";
      if (j.contains("mapvariables")) {
        // std::vector<std::string> vars = j["mapvariables"];
        nlohmann::json vars = j["mapvariables"];
        for (auto &c : vars) {
          if (!c.is_object()) {
            std::string cv = c;
            connection = connection + "\n<map_variables variable_1 = \"" + cv +
                         "\" variable_2 = \"" + cv + "\" />";
          } else { // This is a time variable
            for (auto k : c.items()) {
              timeVariable = k.key();
              break;
            }
          }
        }
      }
      if (j.contains("derivative")) {
        derivativeMapped = true;
        std::string dn = j["derivative"];
        connection = connection + "\n<map_variables variable_1 = \"" + dn +
                     "\" variable_2 = \"dot_" + vname + "\" />";
      }
      // connection = connection + "</connection>";
      return std::make_tuple(filename + ":" + import,
                             compartment + ":" + referenceName + ":" +
                                 connection,
                             "", vname, timeVariable, derivativeMapped);
    } else {
      std::string compartment = j["compartment"];
      std::string targetVar;
      if (j["target"].is_object()) { // User defined
        std::map<std::string, std::string> targetVars = j["target"];
        targetVar = targetVars[vname];
      } else {
        targetVar = j["target"];
      }

      std::string referenceName = j["importName"];
      std::string filename = j["filename"];
      std::replace(filename.begin(), filename.end(), ' ', '_');
      std::string import = "<import xlink:href = \"" + filename + "\" >";
      import = import + "<component component_ref = \"" + compartment +
               "\" name = \"" + referenceName + "_IMP\" /> </import>";
      // std::string connection = "<connection>\n<map_components component_1 =
      // \"" + vname + "_IMP\" component_2 = \"main\" />"; connection =
      // connection + "<map_variables variable_1 = \"" + targetVar + "\"
      // variable_2 = \"" + vname + "\" />";
      std::string connection = "<map_variables variable_1 = \"" + targetVar +
                               "\" variable_2 = \"" + vname + "\" />";
      if (j.contains("stateValue") &&
          j["stateValue"]) { // State variable mapping
        // in this case map to a variable named vname_smap, which will be
        // provided as init for the state variable
        connection = "<map_variables variable_1 = \"" + targetVar +
                     "\" variable_2 = \"" + vname + "_smap\" />";
      }

      std::string timeVariable = "";
      if (j.contains("mapvariables")) {
        nlohmann::json vars = j["mapvariables"];
        for (auto c : vars) {
          if (!c.is_object()) {
            std::string cv = c;
            connection = connection + "\n<map_variables variable_1 = \"" + cv +
                         "\" variable_2 = \"" + cv + "\" />";
          } else { // This is a time variable
            for (auto k : c.items()) {
              timeVariable = k.key();
              break;
            }
          }
        }
      }
      bool derivativeMapped = false;
      if (j.contains("derivative")) {
        derivativeMapped = true;
        std::map<std::string, std::string> dn = j["derivative"];
        for (auto c : dn) {
          connection = connection + "\n<map_variables variable_1 = \"" +
                       c.first + "\" variable_2 = \"" + c.second + "\" />";
        }
      }
      // connection = connection + "</connection>";
      return std::make_tuple(filename + ":" + import,
                             compartment + ":" + referenceName + ":" +
                                 connection,
                             "", vname, timeVariable, derivativeMapped);
    }
  } catch (std::exception &e) {
    throw BGException("Failed to generate cellml linkages for " + vname +
                      " with value " + value + " (A valid JSON is expected)");
  }
}

std::map<std::string, std::string>
getCellML(std::string modelName_, const RCPLIB::RCP<BondGraphInterface> &host_,
          ComputeEquationResults &bgequations) {
  std::string modelName = modelName_;
  std::replace(modelName.begin(), modelName.end(), ' ', '_');
  const RCPLIB::RCP<BondGraph> host =
      RCPLIB::rcp_dynamic_cast<BondGraph>(host_);
  bool solvable = bgequations.bondGraphValidity;
  std::map<RCPLIB::RCP<const SymEngine::Basic>,
           RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
      &equations = bgequations.dof;
  std::map<RCPLIB::RCP<const SymEngine::Basic>,
           RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
      &dofConstrains = bgequations.dof_constraints;
  std::map<RCPLIB::RCP<const SymEngine::Basic>,
           RCPLIB::RCP<const SymEngine::Basic>, SymEngine::Comparator>
      &bondEquations = bgequations.bondEquations;
  SymEngine::vec_basic &constraints = bgequations.constraints;
  std::unordered_map<std::string, std::tuple<std::string, std::string, char>>
      dimensions = bgequations.physicalDimensions;
  std::unordered_map<std::string,std::vector<nlohmann::json>>& annotations = bgequations.annotations;
  std::map<std::string, std::string> files;
  std::vector<std::string> connections;
  std::map<std::string, std::string> imports;
  std::vector<std::string> blocks;
  std::map<std::string, std::vector<unsigned int>> cmetaid;
  unsigned int metaid = 1;

  std::ostringstream cellML, metaidGen,cmetass;
  // std::tie(solvable, equations, bondEquations, constraints, dimensions) =
  // bgequations; Find equations whose lhs has dot_
  std::unordered_map<std::string, std::string> derivatives;
  for (auto eq : equations) {
    cellML.str("");
    cellML.clear();
    cellML << *eq.first;
    std::string variable = cellML.str();
    if (variable.rfind("dot_", 0) == 0) {
      derivatives[variable.substr(4)] = variable;
    }
  }
  cellML.str("");
  cellML.clear();
  std::ostringstream mapParameterFile;
  std::ostringstream connect;

  mapParameterFile << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
  mapParameterFile
      << "<model name=\"" << modelName
      << "Parameters \" xmlns=\"http://www.cellml.org/cellml/1.1#\" "
         "xmlns:cellml=\"http://www.cellml.org/cellml/1.1#\" "
         "xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
      << std::endl;
  CellMLUnits unitMap;
  std::map<std::string, bool> mapUnits;
  std::map<std::string, units::precise_unit> existingUnits;
  std::map<std::string, std::string> defineUnits;
  std::map<std::string, std::string> definedUnitNames;
  std::map<std::string, std::string> dimUnitMap;
  std::map<std::string,
           std::tuple<std::string, std::string, units::precise_unit>>
      cellMLDefinitions;
  // Create variables and their annotations
  std::map<std::string, std::tuple<std::string, std::string, char>> orderedDim(
      dimensions.begin(), dimensions.end());

  auto getMetaIDString = [&metaidGen](unsigned int mi){
    metaidGen.str("");
    metaidGen.clear();
    metaidGen <<"id_"<<std::setw(9) << std::setfill('0')<< mi;
    return metaidGen.str();
  };

  std::map<std::string, std::tuple<std::string, std::string, char>> newDims;
  for (auto var : orderedDim) {
    std::string dim = std::get<0>(var.second);
    std::string unitName = dim;

      if (dim != "") {       
        unitName = unitMap.getUnitName(dim);
        if(unitName!="UNIT_NAME_NOT_FOUND"){
          definedUnitNames[var.first] = unitName;
          mapUnits[unitName] = true;
          dimUnitMap[dim] = unitName;
        }else{
          newDims[var.first] = var.second;
        }
      } else {
        unitName = "dimensionless";
        definedUnitNames[var.first] = unitName;
        // mapUnits[unitName] = true; //dimensionless is predefined, so no need
        // to define it
        existingUnits[unitName] = units::precise::one;
        dimUnitMap[dim] = unitName;
      }
  }
  for (auto var : newDims) {
    std::string dim = std::get<0>(var.second);
    std::string unitName = dim;
    std::tuple<std::string, std::string, units::precise_unit> res;
    if (cellMLDefinitions.find(dim) == cellMLDefinitions.end()) {
      res = unitMap.getCellMLDef(unitName);
      cellMLDefinitions[dim] = res;
    } else {
      res = cellMLDefinitions[dim];
    }
    // auto res = unitMap->getCellMLDef(unitName);
    units::precise_unit unitDef = std::get<2>(res);
    bool found = false;
    // Check if unit definition already exists
    for (auto &b : existingUnits) {
      if (b.second == unitDef) {
        unitName = b.first;
        found = true;
        break;
      }
    }
    if (!found) {
      existingUnits[unitName] = unitDef;
      defineUnits[unitName] = std::get<1>(res);
    }
    definedUnitNames[var.first] = std::get<0>(res);
    dimUnitMap[dim] = unitName;
  }

  // Define cellml units
  // Map known units
  {
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLPrinter printer;
    std::string units = unitMap.getPredefinedunitsXML();
    if (doc.Parse(units.c_str()) == tinyxml2::XML_SUCCESS) {
      doc.Print(&printer);
      files["Units.cellml"] = printer.CStr();
    } else {
      files["Units.cellml"] = units;
      logWarn("Failed to parse standard units definition file!!");
    }
  }

  connect.str("");
  connect.clear();
  connect << "<import xlink:href = \"Units.cellml\">" << std::endl;
  // mapParameterFile << "<import xlink:href = \"Units.cellml\">" << std::endl;
  for (auto mu : mapUnits) {
    connect << "<units name = \"" << mu.first << "\" units_ref=\"" << mu.first
            << "\" />" << std::endl;
    // mapParameterFile << "<units name = \"" << mu.first << "\" units_ref=\""
    // << mu.first << "\" />" << std::endl;
  }
  connect << "</import>" << std::endl;
  // mapParameterFile << "</import>" << std::endl;
  {
    std::string reS = connect.str();
    connections.push_back(reS); // Units
  }
  connect.str("");
  connect.clear();
  // Setup all units in main file
  for (auto un : defineUnits) {
    connect << un.second << std::endl;
  }
  // Setup parameter units in parameter file
  std::map<std::string, bool> definedUnits;
  for (auto var : orderedDim) {
    if (std::get<2>(var.second) == 'p') {
      std::string dim = std::get<0>(var.second);
      // Check if value is json, in which case another map is created
      auto js = std::get<1>(var.second).find('{') == std::string::npos;
      if (definedUnits.find(dim) == definedUnits.end() && js) {
        auto unitName = dimUnitMap[dim];
        mapParameterFile << defineUnits[unitName] << std::endl;
        definedUnits[dim] = true;
      }
    }
  }

  {
    std::string reS = connect.str();
    blocks.push_back(reS); // Units
  }
  connect.str("");
  connect.clear();
  // Define link to parameters file
  connect << "<import xlink:href = \"" << modelName << "Parameters.cellml\" >"
          << std::endl;
  connect
      << "<component component_ref = \"parameters\" name = \"parameters\" />"
      << std::endl;
  connect << "</import>" << std::endl;
  {
    std::string reS = connect.str();
    imports["parameters"] = reS; // Units
  }
  connect.str("");
  connect.clear();
  // Define all the parameters and map them to cellML
  connect << "<connection>\n<map_components component_1 = \"parameters\" "
             "component_2 = \""
          << modelName << "_main\" />" << std::endl;

  mapParameterFile << "<component name=\"parameters\">" << std::endl;
  for (auto var : orderedDim) {
    if (std::get<2>(var.second) == 'p') {
      std::string dim = std::get<0>(var.second);
      // Check if value is json, in which case another map is created
      auto js = std::get<1>(var.second).find('{') == std::string::npos;
      if (js) {
        std::string unitName = definedUnitNames[var.first];
        connect << "<map_variables variable_1 = \"" << var.first
                << "\" variable_2 = \"" << var.first << "\" />" << std::endl;
        mapParameterFile << "<variable initial_value = \""
                         << std::get<1>(var.second) << "\" name = \""
                         << var.first
                         << "\" public_interface = \"out\" units = \""
                         << unitName << "\" />" << std::endl;
      }
    }
  }
  connect << "</connection>" << std::endl;
  mapParameterFile << "</component>\n</model>" << std::endl;
  {
    std::string reS = connect.str();
    connections.push_back(reS); // Units
  }
  {
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLPrinter printer;
    std::string params = mapParameterFile.str();
    if (doc.Parse(params.c_str()) == tinyxml2::XML_SUCCESS) {
      doc.Print(&printer);
      files[modelName + "Parameters.cellml"] = printer.CStr();
    } else {
      files[modelName + "Parameters.cellml"] = params;
      logWarn("Failed to parse generated parameters definition file!!");
    }
  }
  std::unordered_map<std::string, bool> computedVariables;
  std::unordered_map<std::string, std::string> mappedVariables;
  for (auto eq : equations) {
    std::string eqs = eq.first->__str__();
    computedVariables[eqs] = true;
  }
  for (auto eq : bondEquations) {
    std::string eqs = eq.first->__str__();
    computedVariables[eqs] = true;
  }
  // All variables that appear in a constraint are computed i.e. no initial
  // value should be provided
  for (auto c : constraints) {
    SymEngine::set_basic atoms =
        SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*c);
    for (auto eq : atoms) {
      std::string eqs = eq->__str__();
      computedVariables[eqs] = true;
    }
  }
  // Find all variables that need to be mapped
  // Control variables, userdefined variables - values must be JSON
  std::map<std::string, std::string> importMap;
  std::map<std::string, std::vector<std::tuple<std::string, std::string>>>
      connectionMap;
  std::map<std::string, std::string> timeVariableMap;

  for (auto var : orderedDim) {
    std::string value = std::get<1>(var.second);
    // Check if json
    if (value.find('{') != std::string::npos) {
      // Generate import, connection and file for value and the name to import
      auto g = generateLinkages(var.first, var.second);

      std::string imp = std::get<0>(g);
      auto ss = imp.find(":<");
      auto importFile = imp.substr(0, ss);
      // Cut to basename
      //  std::string baseFilename =
      //  importFile.substr(importFile.find_last_of("/\\") + 1);
      //  std::string::size_type const p(baseFilename.find_last_of('.'));
      //  std::string filenameWithoutExtension = baseFilename.substr(0, p);
      //  importMap[filenameWithoutExtension] = imp.substr(ss + 1);
      std::string imf = imp.substr(ss + 1);

      imp = std::get<1>(g);
      ss = imp.find(':');
      auto compartment = imp.substr(0, ss);
      auto nxbit = imp.substr(ss + 1);
      auto si = nxbit.find(':');
      auto variable = nxbit.substr(0, si);
      importMap[variable] = imf;
      // Time variable map
      if (std::get<4>(g) != "") {
        timeVariableMap[variable + "_IMP"] = std::get<4>(g);
      }
      if (connectionMap.find(compartment) != connectionMap.end()) {
        connectionMap[compartment].push_back(
            std::make_tuple(variable, nxbit.substr(si + 1)));
      } else {
        connectionMap[compartment] = {
            std::make_tuple(variable, nxbit.substr(si + 1))};
      }
      if (std::get<2>(g) != "")
        files[var.first + ".cellml"] = std::get<2>(g);
      mappedVariables[var.first] = std::get<3>(g);
      // Add the name to computedVariables to not generate the initialValue
      // computedVariables[var.first] = true;
      if (std::get<2>(var.second) == 'c' &&
          std::get<5>(g)) { // Create a map only if a mapping is available
        mappedVariables["dot_" + var.first] = "dot_" + std::get<3>(g);
      }
    }
  }
  for (auto c : importMap) {
    imports[c.first + "_IMP"] = c.second;
  }
  for (auto c : timeVariableMap) {
    std::string connection = "<connection>\n<map_components component_1 = \"" +
                             c.first + " \" component_2 =\"time\" />";
    connection = connection + "<map_variables variable_1=\"" + c.second +
                 "\" variable_2=\"t\"/>";
    connection = connection + "</connection>";
    connections.push_back(connection);
  }
  for (auto c : connectionMap) {
    std::map<std::string, std::vector<std::string>> vmap;
    for (auto v : c.second) {
      if (vmap.find(std::get<0>(v)) != vmap.end()) {
        vmap[std::get<0>(v)].push_back(std::get<1>(v));
      } else {
        vmap[std::get<0>(v)] = {std::get<1>(v)};
      }
    }
    std::string comp = c.first;
    for (auto v : vmap) {
      std::string connection =
          "<connection>\n<map_components component_1 = \"" + v.first +
          "_IMP\" component_2 =\"" + modelName + "_main\" />";
      for (auto x : v.second) {
        connection = connection + x;
      }
      connection = connection + "</connection>";
      connections.push_back(connection);
    }
  }
  // Start generation
  cellML.str("");
  cellML.clear();
  cellML << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
  cellML << "<model name=\"" << modelName
         << "\" xmlns=\"http://www.cellml.org/cellml/1.1#\" "
            "xmlns:cellml=\"http://www.cellml.org/cellml/1.1#\" "
            "xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
         << std::endl;
  // Define units
  cellML << connections[0] << std::endl;
  connections.erase(connections.begin()); // Remove units import
  cellML << blocks[0] << std::endl;
  blocks.erase(blocks.begin()); // Remove units
  // Define imports
  std::ostringstream encapsulationDef;

  encapsulationDef
      << "<group><relationship_ref "
         "relationship=\"encapsulation\"/><component_ref component=\""
      << modelName << "_main\">";
  // Not adding time component to encapsulation to break assignment
  // encapsulationDef << "<component_ref component=\""<<modelName<<"_time\"/>";
  for (auto c : imports) {
    encapsulationDef << "<component_ref component=\"" << c.first << "\"/>";
    cellML << c.second << std::endl;
  }
  encapsulationDef << "</component_ref></group>";
  // Define connections
  for (auto c : connections) {
    cellML << c << std::endl;
  }
  // Time connection
  cellML << "<connection><map_components component_1=\"time\" component_2=\""
         << modelName
         << "_main\"/><map_variables variable_1=\"t\" "
            "variable_2=\"t\"/></connection>";
  // Define other blocks
  for (auto c : blocks) {
    cellML << c << std::endl;
  }
  // Define encapsulation
  cellML << encapsulationDef.str() << std::endl;
  // Define time component
  cellML << "<component name=\"time\">  <variable name=\"t\" "
            "public_interface=\"out\" units=\"second\"/> </component>";
  // Define the component
  cellML << "<component name=\"" << modelName << "_main\">" << std::endl;

  if (!solvable) {
    cellML
        << "<!-- The bondgraph library has identified issues with this model. "
           "The cellml description may not build or give correct results, view "
           "the library's log (if enabled) for more details -->"
        << std::endl;
  }
  // Define time variable
  cellML << "<variable name=\"t\" units=\"second\" public_interface = \"in\"/>"
         << std::endl;

  for (auto var : orderedDim) {
    std::string dim = std::get<0>(var.second);
    std::string unitName = definedUnitNames[var.first];
    std::string value = std::get<1>(var.second);
    std::string initVal = "initial_value = \"" + value + "\"";
    if (computedVariables[var.first]) {
      initVal = "";
    }


    //Check here if a variable has an annotation and add to the variable definition
    //<variable cmeta:id="id_00000000" initial_value="" name="" units=""/>
    
    cmetass.str("");
    cmetass.clear();
    if(annotations.find(var.first)!=annotations.end()){
       for(size_t i=0;i<annotations[var.first].size();i++){
        if(i==0){
          cmetass<<getMetaIDString(metaid);
        }else{
          cmetass<<" id=\""<<getMetaIDString(metaid)<<"\"";
        }
        cmetaid[var.first].push_back(metaid++);
       }
    }
    std::string cmeta = cmetass.str();
    switch (std::get<2>(var.second)) {
    case 'p': {
      if(cmeta.size()==0){
        cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName
              << "\" public_interface = \"in\" />";
      }else{
        cellML << "< variable cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\"" << unitName
              << "\" public_interface = \"in\" xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";
      }
      break;
    }
    case 'c': {
      // TODO: These need to be mapped to a component
      if (mappedVariables.find(var.first) == mappedVariables.end()) {
        if(cmeta.size()==0){
          cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName
               << "\" " << initVal << "/>";
        }else{
          cellML << "< variable cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\"" << unitName
               << "\" " << initVal << " xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";
        }
      } else {
        if(cmeta.size()==0){
          cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName
                << "\" public_interface = \"in\" />";
        }else{
          cellML << "< variable cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\"" << unitName
                << "\" public_interface = \"in\" xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";
        }
      }
      break;
    }
    case 's': {
      if (var.first.size() > 4 &&
          var.first.substr(0, 4) == "dot_") { // Derivatives are computed
        if(cmeta.size()==0){
          cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName
                << "\" public_interface = \"out\" />";
        }else{
          cellML << "< variable cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\"" << unitName
                << "\" public_interface = \"out\" xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";
        }
      } else {
        if (value.find(":") == std::string::npos) // If value is not json
          if(cmeta.size()==0){
            cellML << "< variable  name=\"" << var.first << "\" units=\""
                  << unitName << "\" public_interface = \"out\" " << initVal
                  << "/>";
          }else{
            cellML << "< variable  cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\""
                  << unitName << "\" public_interface = \"out\" " << initVal
                  << "xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";            
          }
        else { // The value is json and is mapped
          // Create a variable vname_smap that will be mapped
          // Set that variable vname_smap as the initial value
          cellML << "< variable  name=\"" << var.first << "_smap\" units=\""
                 << unitName << "\" public_interface = \"in\" />";
          if(cmeta.size()==0){
            cellML << "< variable  name=\"" << var.first << "\" units=\""
                  << unitName << "\" public_interface = \"out\""
                  << "initial_value = \"" + var.first + "_smap\" />";
          }else{
            cellML << "< variable  cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\""
                  << unitName << "\" public_interface = \"out\""
                  << "initial_value = \"" + var.first + "_smap\" xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";
          }
        }
      }
      break;
    }
    default: { // Variables
      if(cmeta.size()==0){
        cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName
              << "\" public_interface = \"out\" />";
      }else{
        cellML << "< variable cmeta:id=\""<<cmeta<<"\" name=\"" << var.first << "\" units=\"" << unitName
              << "\" public_interface = \"out\" xmlns:cmeta=\"http://www.cellml.org/metadata/1.0#\" />";
      }
      break;
    }
    }
  }
  cellML << "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">" << std::endl;

  // Create Bond variable equations
  std::vector<std::string> bondEqMathML = getMathML(bondEquations);
  cellML << "<!--Port values -->" << std::endl;
  for (auto expr : bondEqMathML) {
    cellML << expr << std::endl;
  }
  // Create Constraint equations
  if (constraints.size() > 0) {
    cellML << "<!--Constraint equations -->" << std::endl;
    std::vector<std::string> consEqMathML = getMathML(constraints);
    for (auto expr : consEqMathML) {
      cellML << expr << std::endl;
    }
  }

  std::vector<std::string> stateEqMathML = getMathML(equations);
  cellML << "<!--State equations -->" << std::endl;
  for (auto expr : stateEqMathML) {
    cellML << expr << std::endl;
  }
  if (dofConstrains.size() > 0) {
    cellML << "<!-- List of state variables that have an algebriac solution, "
              "select the ode or algebraic version as desired -->"
           << std::endl;
    std::vector<std::string> stateSolEqMathML = getMathML(dofConstrains);
    for (auto expr : stateSolEqMathML) {
      cellML << "<!-- " << expr << " -->" << std::endl;
    }
    // std::map<std::string, bool> skipEq;
    // for (auto dc : dofConstrains) {
    //     skipEq[dc.first->__str__()] = true;
    // }
    // for (auto expr : derivatives) {
    //     // Provide ode if state equation solution is not available
    //     if (skipEq.find(expr.first) == skipEq.end())
    //         cellML << "<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>"
    //         << expr.first << "</ci></apply><ci>" << expr.second <<
    //         "</ci></apply> " << std::endl;
    //     else
    //         cellML << "<!--"
    //                << "<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>"
    //                << expr.first << "</ci></apply><ci>" << expr.second <<
    //                "</ci></apply> "
    //                << "-->" << std::endl;
    // }
    for (auto expr : derivatives) {
      cellML << "<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>"
             << expr.first << "</ci></apply><ci>" << expr.second
             << "</ci></apply> " << std::endl;
    }
    std::ostringstream solx;
    solx << dofConstrains.size();
    files["StateEquationSolutions"] = solx.str();
  } else {
    for (auto expr : derivatives) {
      cellML << "<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>"
             << expr.first << "</ci></apply><ci>" << expr.second
             << "</ci></apply> " << std::endl;
    }
  }
  cellML << "</math>" << std::endl;
  cellML << "</component>" << std::endl;
  //Generate RDF for annotations
  if(cmetaid.size()>0){
    cellML << "<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\">" << std::endl;
    for(const auto& mi : cmetaid){
      std::vector<nlohmann::json>& avec = annotations[mi.first];   
      const std::vector<unsigned int>& ids = mi.second; 
      cellML <<"\t<!--"<<mi.first<<"-->"<<std::endl;
      for(size_t i=0;i<avec.size();i++){
        auto& anot = avec[i];
        std::string metaid = getMetaIDString(ids[i]);
        cellML << "\t<rdf:Description rdf:about=\"#"<<metaid<<"\">"<<std::endl;
        std::string rel = anot["relationship"];
        std::string identifiers_org_uri = anot["annotation"]["identifiers_org_uri"];
        auto ploc = rel.rfind(":");
        if (ploc != std::string::npos) {
          rel = rel.substr(ploc+1,rel.size());
        }         
        cellML << "\t\t<"<<rel<<" xmlns=\"http://biomodels.net/biology-qualifiers/\">"<<std::endl;
        cellML << "\t\t\t<rdf:Description rdf:about=\""<<identifiers_org_uri<<"\"/>"<<std::endl;
        cellML << "\t\t</"<<rel<<">"<<std::endl;
        cellML << "\t</rdf:Description>"<<std::endl;
      }
    }
    cellML << "</rdf:RDF>" << std::endl;
  }
  cellML << "</model>";

  std::string computedCellML = cellML.str();
  // Do name mapping

  tinyxml2::XMLDocument doc;
  tinyxml2::XMLPrinter printer;

  if (doc.Parse(computedCellML.c_str()) == tinyxml2::XML_SUCCESS) {
    doc.Print(&printer);
    files[modelName + ".cellml"] = printer.CStr();
    return files;
  } else {
    logDebug(computedCellML);
    // std::cout<<"\n\n"<<computedCellML<<"\n\n"<<std::endl;
    throw BGException("Error creating CellML\n" + cellML.str());
  }
}
} // namespace BG