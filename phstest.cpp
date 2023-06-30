#include "src/bondgraph.hpp"
#include "src/thirdparty/json.hpp"
#include <fstream>

using namespace BG;

RCPLIB::RCP<BondGraphInterface> loadProject(std::string file,
                                            bool phs = false) {
  std::ifstream ifs(file);
  nlohmann::json jf = nlohmann::json::parse(ifs);
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
      // std::replace(displayName.begin(), displayName.end(), ':', 'c');
      std::replace_if(
          displayName.begin(), displayName.end(),
          [](auto ch) { return std::ispunct(ch); }, '_');
      std::replace(displayName.begin(), displayName.end(), ' ', '_');

      auto cellmlName = displayName;
      auto dName = def["mid"]; // Connections store mid's

      std::string domain = edef["domain"];
      if (domain == "Annotation") {
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

  std::cout << file << std::endl;
  if (phs) {
    ioBondGraph->computeStateEquationNoDim();
    // nlohmann::json res = ioBondGraph->computePortHamiltonian();
    // std::cout << res.dump() << std::endl;
  } else {
    auto eqs = ioBondGraph->computeStateEquation();
    auto files = getCellML("RLC", ioBondGraph, eqs);
    std::cout << files["RLC.cellml"] << std::endl;
  }
  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> simpleRCV() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createCapacitor();
  lC1->setParameter("C", "1.0", "mega Farad");
  ioBondGraph->addComponent(lC1);

  // Create the resistor
  auto lR = createResistor();
  ioBondGraph->addComponent(lR);

  // Create the junctions
  auto lJ1_1 = createOneJunction();
  ioBondGraph->addComponent(lJ1_1);

  // Create Source
  auto lSe = createConstantVoltageSource();
  ioBondGraph->addComponent(lSe);

  // Create the bonds
  ioBondGraph->connect(lJ1_1, lR);
  ioBondGraph->connect(lJ1_1, lC1);
  ioBondGraph->connect(lJ1_1, lSe);

  std::cout << " Simple RC+V " << std::endl;
  nlohmann::json res = ioBondGraph->computePortHamiltonian();
  std::cout << res.dump() << std::endl;
  auto eqs = ioBondGraph->computeStateEquation();
  auto files = getCellML("RCV", ioBondGraph, eqs);
  std::cout << files["RCV.cellml"] << std::endl;
  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> rlc() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createCapacitor();
  ioBondGraph->addComponent(lC1);

  auto lC2 = createCapacitor();
  ioBondGraph->addComponent(lC2);

  // Create the Flow source
  auto lSf = createConstantFlowSource();
  ioBondGraph->addComponent(lSf);

  // Create the resistor
  auto lR = createResistor();
  ioBondGraph->addComponent(lR);
  auto lR2 = createResistor();
  ioBondGraph->addComponent(lR2);

  // Create the Transformer
  auto lTf = createTransformer();
  ioBondGraph->addComponent(lTf);

  // Create the junctions
  auto lJ0_1 = createZeroJunction();
  auto lJ1_1 = createOneJunction();
  ioBondGraph->addComponent(lJ0_1);
  ioBondGraph->addComponent(lJ1_1);

  // Create the bonds
  ioBondGraph->connect(lJ1_1, lR);
  ioBondGraph->connect(lJ1_1, lC1);
  ioBondGraph->connect(lJ1_1, lTf);

  ioBondGraph->connect(lTf, 1, lJ0_1);
  ioBondGraph->connect(lJ0_1, lC2);
  ioBondGraph->connect(lJ0_1, lR2);
  ioBondGraph->connect(lSf, lJ0_1);

  std::cout << " PHS " << std::endl;
  nlohmann::json res = ioBondGraph->computePortHamiltonian();
  std::cout << res.dump() << std::endl;
  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> reaction() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto A = createConcentration();
  A->setName("A");
  ioBondGraph->addComponent(A);

  auto B = createConcentration();
  B->setName("B");
  ioBondGraph->addComponent(B);

  // Create the reaction
  auto Re = createReaction();
  Re->setName("Re");
  ioBondGraph->addComponent(Re);

  // Create the junctions
  auto Y_A = createOneJunction();
  Y_A->setName("Y_A_1");
  ioBondGraph->addComponent(Y_A);

  auto Y_B = createOneJunction();
  Y_B->setName("Y_B_1");
  ioBondGraph->addComponent(Y_B);

  // Create the bonds
  ioBondGraph->connect(A, Y_A);
  ioBondGraph->connect(B, Y_B);
  ioBondGraph->connectInverting(Re, 0, Y_A);
  ioBondGraph->connectInverting(Re, 1, Y_B);

  std::cout << " Reaction " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> eReaction() {
  // Model reaction using electrical elements
  // Reaction is modelled as 1 junction with a resistor
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createCapacitor();
  lC1->setName("A");
  ioBondGraph->addComponent(lC1);

  auto lC2 = createCapacitor();
  lC2->setName("B");
  ioBondGraph->addComponent(lC2);

  // Create the resistor
  auto lR = createResistor();
  lR->setName("R");
  ioBondGraph->addComponent(lR);

  // Create the junctions
  auto lJ1_A = createOneJunction();
  auto lJ1_B = createOneJunction();
  auto lJ1_Re = createOneJunction();

  ioBondGraph->addComponent(lJ1_A);
  ioBondGraph->addComponent(lJ1_B);
  ioBondGraph->addComponent(lJ1_Re);

  // Create the bonds
  ioBondGraph->connect(lJ1_Re, lR);
  ioBondGraph->connect(lJ1_Re, lJ1_A);
  ioBondGraph->connect(lJ1_Re, lJ1_B);
  ioBondGraph->connect(lJ1_A, lC1);
  ioBondGraph->connect(lJ1_B, lC2);

  std::cout << " eReaction " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

int main(int argc, char *argv[]) {
  // rlc();
   simpleRCV();
   reaction();
   eReaction();
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/GPCRC/GPCRReactionC.json");
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Demonstration/Demonstration.json");
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/RC/RCcircuitWUI.json");
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Test/bve.json");
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Test/Fail1.json");
  //loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Test/DimCheck.json");
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Test/Memristor.json",true);
  // loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Test/Composite.json",true);
  loadProject("/mnt/d/GithubRepositories/BGUITest/Examples/Test/RC2Units.json");//,
  //             true);
  return 0;
}