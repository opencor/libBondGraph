//#define DEBUG_PHS
#define CHECK_DEPENDEND_SOURCES
#include <fstream>

static int getMatrixRank(SymbolicMatrix &mat) {
  Eigen::MatrixX<float> R(mat.rows(), mat.cols());
  R.fill(0.0);
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      if (!(mat(i, j).get_basic() == SymEngine::zero))
        R(i, j) = 1.0;
    }
  }
  Eigen::FullPivLU<Eigen::MatrixX<float>> lu_decomp(R);
  return lu_decomp.rank();
};

// Function gives B (i) < - and B (i) ->
static std::tuple<
    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>,
    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>>
giveBondsInOut(
    RCPLIB::RCP<BGElement> &element,
    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
        &setBond,
    std::map<std::string, unsigned int> &elementIndexs) {
  std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
      bondIn;
  std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
      bondOut;

  for (auto j = 0; j < setBond.size(); j++) {
    RCPLIB::RCP<BGElement> &start = std::get<0>(setBond.at(j));
    RCPLIB::RCP<BGElement> &finish = std::get<1>(setBond.at(j));
    auto eid = elementIndexs[element->getId()];
    auto sid = elementIndexs[start->getId()];
    auto fid = elementIndexs[finish->getId()];
#ifdef DEBUG_PHS
    std::cout << "Considering " << eid << " BGin " << sid << "x" << fid
              << std::endl;
#endif
    if (start->getId() == element->getId()) { // Bond is ingoing in element
      bondOut.push_back(setBond.at(j));
#ifdef DEBUG_PHS
      std::cout << eid << " BGin " << sid << "x" << fid << std::endl;
#endif
    } else if (finish->getId() ==
               element->getId()) { // Bond is outgoing in element
      bondIn.push_back(setBond.at(j));
#ifdef DEBUG_PHS
      std::cout << eid << " BGout " << sid << "x" << fid << std::endl;
#endif
    }
  }
  return std::make_tuple(bondIn, bondOut);
};

// Function gives set of incident bonds B(i)
static std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
giveBondsIncident(
    RCPLIB::RCP<BGElement> &element,
    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
        &setBond) {
  std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>> setBi;
  for (auto j = 0; j < setBond.size(); j++) {
    RCPLIB::RCP<BGElement> &start = std::get<0>(setBond.at(j));
    RCPLIB::RCP<BGElement> &finish = std::get<1>(setBond.at(j));
    if (start->getId() == element->getId() ||
        finish->getId() ==
            element->getId()) { // bond j is incident to element i
      setBi.push_back(setBond.at(j));
    }
  }
  return setBi;
};

static int getEigenVectorIndex(SymEngine::Expression &elem,
                               SymbolicMatrix &vecFIC) {
  for (int x = 0; x < vecFIC.rows(); x++) {
    if (vecFIC(x, 0) == elem) {
      return x;
    }
  }
  return -1;
};

nlohmann::json to_json(SymbolicMatrix &mat) {
  nlohmann::json js;
  js["rows"] = mat.rows();
  js["cols"] = mat.cols();
  std::vector<std::string> elems;
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      SymEngine::Expression simpE = SymEngine::simplify(mat(i, j));
      // elems.push_back(SymEngine::mathml(simpE));
      elems.push_back(SymEngine::latex(simpE));
    }
  }
  js["elements"] = elems;
  return js;
}

nlohmann::json BondGraph::computePortHamiltonian() {
  nlohmann::json result; // Compose results as we progress
  result["success"] = false;
  std::vector<std::string> warnings;
  try {
    std::vector<RCPLIB::RCP<BGElement>> connectedComponents;
    std::map<std::string, RCPLIB::RCP<BGElement>> compIDMap;
    const SymEngine::Expression symNegOne = SymEngine::Expression("-1");
    const SymEngine::Expression symOne = SymEngine::Expression("1");
    const SymEngine::Expression symZero = SymEngine::Expression("0");

    for (auto &bd : mBonds) {
      auto from = bd->getFromPort()->getComponent();
      auto to = bd->getToPort()->getComponent();
      compIDMap[from->getId()] = from;
      compIDMap[to->getId()] = to;
    }

    // Collect different types of elements
    std::vector<RCPLIB::RCP<BGElement>> setVC;     // Storages
    std::vector<RCPLIB::RCP<BGElement>> setVConc;  // Concentrations
    std::vector<RCPLIB::RCP<BGElement>> setVR;     // Dissipative elements
    std::vector<RCPLIB::RCP<BGElement>> setVRe;    // Reactive elements
    std::vector<RCPLIB::RCP<BGElement>> setVSf;    // Flow sources
    std::vector<RCPLIB::RCP<BGElement>> setVSe;    // Effort sources
    std::vector<RCPLIB::RCP<BGElement>> setV0;     // 0-Juntions
    std::vector<RCPLIB::RCP<BGElement>> setV1;     // 1-Juntions
    std::vector<RCPLIB::RCP<BGElement>> setVTF;    // Transformers
    std::vector<RCPLIB::RCP<BGElement>> setVStoic; // Stoichiometry
    std::vector<RCPLIB::RCP<BGElement>> setVGY;    // Gyrators

    for (auto &c : mComponents) { // Maintain order
      auto id = c->getId();

      if (compIDMap.find(id) != compIDMap.end()) {
        connectedComponents.push_back(compIDMap[id]);
        auto &mc = compIDMap[id];
        switch (mc->getType()) {
        case eCapacitance:
        case eInductance: {
          setVC.push_back(mc);
          break;
        }
        case bConcentration: {
          setVConc.push_back(mc);
          break;
        }
        case bReaction: {
          setVRe.push_back(mc);
          break;
        }
        case eResistance: {
          setVR.push_back(mc);
          break;
        }
        case bStoichiometry: {
          setVStoic.push_back(mc);
          break;
        }
        case eTransformer: {
          setVTF.push_back(mc);
          break;
        }
        case eGyrator: {
          setVGY.push_back(mc);
          break;
        }
        case ePotentialSource: {
          setVSe.push_back(mc);
          break;
        }
        case eFlowSource: {
          setVSf.push_back(mc);
          break;
        }
        case eZero: {
          setV0.push_back(mc);
          break;
        }
        case eOne: {
          setV1.push_back(mc);
          break;
        }
        default: {
          throw BGException(
              "Component of type no supported. Offending component " +
              mc->getName());
        }
        }
      }
    }

    // For reactions move to 1 port representation ( A one junction with the
    // concentrations + a resistor with the reaction rate)
    // See Network thermodynamics of biological systems: A bond graph approach
    // Peter J.Gawthrop Michael Pan
    // Also, Network thermodynamics: dynamic modelling of biophysical systems
    // George F. Oster, Alan S. Perelson and Aharon Katchalsky

    //[setV0, setV1, setVRe (as 1 junction), setVStoic, setVTF, setVGY]; //
    // Interior elements
    std::vector<RCPLIB::RCP<BGElement>> setVI;
    setVI.insert(std::end(setVI), std::begin(setV0), std::end(setV0));
    setVI.insert(std::end(setVI), std::begin(setV1), std::end(setV1));
    setVI.insert(std::end(setVI), std::begin(setVRe), std::end(setVRe));
    setVI.insert(std::end(setVI), std::begin(setVStoic), std::end(setVStoic));
    setVI.insert(std::end(setVI), std::begin(setVTF), std::end(setVTF));
    setVI.insert(std::end(setVI), std::begin(setVGY), std::end(setVGY));

    //[setVC, setVConc, setVR, setVRe (as Resistance), setVSf, setVSe]; //
    // Exterior elements
    std::vector<RCPLIB::RCP<BGElement>> setVE;
    setVE.insert(std::end(setVE), std::begin(setVC), std::end(setVC));
    setVE.insert(std::end(setVE), std::begin(setVConc), std::end(setVConc));
    setVE.insert(std::end(setVE), std::begin(setVR), std::end(setVR));
    setVE.insert(std::end(setVE), std::begin(setVRe), std::end(setVRe));
    setVE.insert(std::end(setVE), std::begin(setVSf), std::end(setVSf));
    setVE.insert(std::end(setVE), std::begin(setVSe), std::end(setVSe));

    //[setVI, setVE]; // All elements
    std::vector<RCPLIB::RCP<BGElement>> setV;
    setV.insert(std::end(setV), std::begin(setVI), std::end(setVI));
    setV.insert(std::end(setV), std::begin(setVE), std::end(setVE));

    // Element count
    const int numNC = setVC.size() + setVConc.size(); // Number of storages
    const int numNR =
        setVR.size() + setVRe.size(); // Number of dissipative elements
    const int numNSf = setVSf.size(); // Number of flow sources
    const int numNSe = setVSe.size(); // Number of effort sources
    const int numNS = numNSf + numNSe;
    const int numNI = setVI.size(); // Number of interior elements
    const int numNE = setVE.size(); // Number of exterior elements
    const int numN = numNI + numNE; // Total number of bond graph elements

    if (numN == 0) {
      throw BGException("Bondgraph must have one or more connected elements");
    }

    std::map<std::string, unsigned int> elementIndexs;
    std::map<std::string, unsigned int> reactionIndexs;
    std::map<std::string, unsigned int> bondIndexs;
    unsigned int ix = 0;
    // Create element ids in the order of element types [0,1,TF,GY,C,R,S]
    auto psets = {setV0, setV1, setVRe, setVTF, setVGY, setVC, setVConc, setVR};
    for (auto &set : psets) {
      for (auto &elem : set) {
        if (elementIndexs.find(elem->getId()) == elementIndexs.end()) {
#ifdef DEBUG_PHS
          std::cout << elem->getId() << "\t" << elem->getName() << "\t" << ix
                    << std::endl;
#endif
          elementIndexs[elem->getId()] = ix++;
        } else {
          // Report error
          logDebug("Found index for element, multiple includes in setV0!!");
        }
      }
    }
    // For VRe, resistors add indexes to reaction indexes map
    for (auto &elem : setVRe) {
      if (reactionIndexs.find(elem->getId()) == reactionIndexs.end()) {
#ifdef DEBUG_PHS
        std::cout << "R " << elem->getId() << "\t" << elem->getName() << "\t"
                  << ix << std::endl;
#endif
        reactionIndexs[elem->getId()] = ix++;
      } else {
        // Report error
        logDebug("Found index for element, multiple includes in setVRe!!");
      }
    }
    auto ssets = {setVSf, setVSe};
    for (auto &set : ssets) {
      for (auto &elem : set) {
        if (elementIndexs.find(elem->getId()) == elementIndexs.end()) {
#ifdef DEBUG_PHS
          std::cout << elem->getId() << "\t" << elem->getName() << "\t" << ix
                    << std::endl;
#endif
          elementIndexs[elem->getId()] = ix++;
        } else {
          // Report error
          logDebug("Found index for element, multiple includes in setV0!!");
        }
      }
    }

    // Store element index mapping
    std::vector<nlohmann::json> elementInfo;
    for (auto &elem : setV) {
      nlohmann::json es;
      if (elementIndexs.find(elem->getId()) != elementIndexs.end()) {
        es["id"] = elem->getId();
        es["suffix"] = elementIndexs[elem->getId()];
        es["name"] = elem->getName();
        elementInfo.push_back(es);
      }
    }

    result["Elements"] = elementInfo;

    SymbolicMatrix stateVariables(setVC.size() + setVConc.size(), 1);
    std::vector<std::string> hamiltonianString;
    int ctr = 0;
    for (auto &ie : setVC) {
      int numStates = ie->getNumStates();
      auto values = ie->values();
      std::string state = std::get<1>(values[numStates - 1])->name +
                          std::to_string(elementIndexs[ie->getId()]);
      stateVariables(ctr++, 0) = SymEngine::parse(state);
      std::string parm = std::get<1>(values[numStates])->name + "_" +
                         std::to_string(elementIndexs[ie->getId()]);
      std::string energy = "(1/2)*(1/" + parm + ")*" + state + "**2";
      hamiltonianString.push_back(energy);
    }
    for (auto &ie : setVConc) {
      int numStates = ie->getNumStates();
      auto values = ie->values();
      std::string state = std::get<1>(values[numStates - 1])->name +
                          std::to_string(elementIndexs[ie->getId()]);
      stateVariables(ctr++, 0) = SymEngine::parse(state);
      std::string parm = std::get<1>(values[numStates + 2])->name + "_" +
                         std::to_string(elementIndexs[ie->getId()]);
      std::string pR, pT;
      auto Rv = std::get<1>(values[numStates]);
      if (Rv->universalConstant) {
        pR = Rv->name;
      } else {
        pR = Rv->name + "_" + std::to_string(elementIndexs[ie->getId()]);
      }
      auto Rt = std::get<1>(values[numStates + 1]);
      if (Rt->universalConstant) {
        pT = Rt->name;
      } else {
        pT = Rt->name + "_" + std::to_string(elementIndexs[ie->getId()]);
      }
      // RTq_0(ln(kq_0)−1)
      std::string energy =
          pR + "*" + pT + "*" + state + "(log(" + parm + "*" + state + ")-1.0)";
      hamiltonianString.push_back(energy);
    }
    result["stateVector"] = to_json(stateVariables);
    std::ostringstream hms;
    for (auto &e : hamiltonianString) {
      hms << e + " + ";
    }
    std::string hmString = hms.str();
    number hmexp = SymEngine::simplify(SymEngine::expand(
        SymEngine::parse(hmString.substr(0, hmString.size() - 3))));

    result["hamiltonianLatex"] = SymEngine::latex(hmexp);
    std::string hs = hmString.substr(0, hmString.size() - 3);
    hs.erase(std::remove(hs.begin(), hs.end(), '*'), hs.end());
    result["hamiltonian"] = hs;
    // Create the adjacency matrix
    Eigen::MatrixXi adjacency = Eigen::MatrixXi(numN, numN);
    adjacency.fill(0);

    for (auto &bd : mBonds) {
      auto from = bd->getFromPort()->getComponent();
      auto to = bd->getToPort()->getComponent();
      auto fidx = elementIndexs[from->getId()];
      auto tidx = elementIndexs[to->getId()];
#ifdef DEBUG_PHS
      std::cout << fidx << "-> " << tidx << std::endl;
#endif
      // Standard bond graph literature [W. Borutzky, Bond Graph Methodology:
      // Development and Analysis of Multidisciplinary Dynamic System Models,
      // Springer, London, 2010, p. 59] in which bonds are incoming to
      // storages and resistors and outgoing from sources of flow and effort
      if (from->getType() == eCapacitance || from->getType() == eInductance ||
          from->getType() == eResistance || from->getType() == bConcentration) {
        // Report that the bond is not consistent
        logDebug("Bond direction is incorrect ");
        auto t = fidx;
        fidx = tidx;
        tidx = t;
      }
      if (to->getType() == ePotentialSource || to->getType() == eFlowSource ||
          to->getType() == bChemostat || to->getType() == bFlowstat) {
        // Report that the bond is not consistent
        logDebug("Bond direction is incorrect ");
        auto t = fidx;
        fidx = tidx;
        tidx = t;
      }
#ifdef DEBUG_PHS
      std::cout << " Bond " << from->getId() << " -> " << to->getId() << " ("
                << fidx << " x " << tidx << " ) " << std::endl;
#endif
      adjacency(fidx, tidx) = 1;
    }
    // Define links for the VRe resistances, note this will be from the VRe 1
    // junction to VRe resistor
    for (auto &c : setVRe) {
      auto fidx = elementIndexs[c->getId()];
      auto tidx = reactionIndexs[c->getId()];
      adjacency(fidx, tidx) = 1;
    }

#ifdef DEBUG_PHS
    std::cout << " Adjacency " << std::endl;
    std::cout << adjacency << std::endl;
#endif
    // Internal Bonds*
    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
        setBI;
    // Part of adjacency matrix covering internal elements
    for (auto i = 0; i < numNI; i++) {
      for (auto j = 0; j < numNI; j++) {
        if (adjacency(i, j) == 1) {
#ifdef DEBUG_PHS
          std::cout << " In " << i << " " << j << std::endl;
#endif
          setBI.push_back(std::make_tuple(setVI.at(i), setVI.at(j)));
        }
      }
    }

    // External Bonds
    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
        setBE;
    // Part of adjacency matrix covering resistive elements and storages

    for (auto j = 0; j < numNC + numNR; j++) {
      for (auto i = 0; i < numNI; i++) {
        if (adjacency(i, numNI + j) == 1) {
#ifdef DEBUG_PHS
          std::cout << " out " << i << " " << numNI + j << std::endl;
#endif
          setBE.push_back(std::make_tuple(setVI.at(i), setVE.at(j)));
        }
      }
    }

    // Part of adjacency matrix covering sources
    for (auto i = 0; i < numNS; i++) {
      for (auto j = 0; j < numNI; j++) {
        auto fx = numN - numNS + i;
        if (adjacency(fx, j) == 1) {
#ifdef DEBUG_PHS
          std::cout << " sources " << i << " " << j << " fx " << fx
                    << std::endl;
          std::cout << setVE.at(numNC + numNR + i)->getName() << " -> "
                    << setVI.at(j)->getName() << std::endl;
#endif
          setBE.push_back(
              std::make_tuple(setVE.at(numNC + numNR + i), setVI.at(j)));
        }
      }
    }

    std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
        setB;

    setB.insert(std::end(setB), std::begin(setBE), std::end(setBE));
    setB.insert(std::end(setB), std::begin(setBI), std::end(setBI));

    unsigned int numMI = setBI.size(); // Number of interior multi bonds m_I
    unsigned int numME = setVE.size(); // Number of exterior multi bonds m_E
    // Definition of Dirac structures for each interior element
    std::map<std::string, std::vector<int>> setDirac_; // Store the flow indexes
    std::map<std::string, std::map<int, SymbolicMatrix>> setDirac;

    unsigned int bidx = 0;
    for (auto &b : setBI) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      if (bondIndexs.find(bnd) == bondIndexs.end()) {
        bondIndexs[bnd] = bidx++;
      }
    }

    for (auto &b : setBE) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance

        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      if (bondIndexs.find(bnd) == bondIndexs.end()) {
        bondIndexs[bnd] = bidx++;
      }
    }

    // for all interior elements
    for (auto &ie : setVI) {
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bondsIn;
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bondsOut;
      std::tie(bondsIn, bondsOut) = giveBondsInOut(ie, setB, elementIndexs);
      std::vector<std::string> flows;
      std::vector<int> flowNum; // The bond number
      std::vector<std::string> efforts;
      std::map<int, SymbolicMatrix> dirac;

      for (auto &b : bondsIn) {
        RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
        RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
        std::string bnd;
        if (e1->getId() != e2->getId()) { // Handle VRe elements which have both
                                          // a 1 junction and a resistance

          bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
                std::to_string(elementIndexs[e2->getId()]);
        } else {
          bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
                std::to_string(reactionIndexs[e2->getId()]);
        }

        flowNum.push_back(bondIndexs[bnd]);
        flows.push_back("f_" + bnd);
        efforts.push_back("e_" + bnd);
      }

      for (auto &b : bondsOut) {
        RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
        RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
        std::string bnd;
        if (e1->getId() != e2->getId()) { // Handle VRe elements which have both
                                          // a 1 junction and a resistance

          bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
                std::to_string(elementIndexs[e2->getId()]);
        } else {
          bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
                std::to_string(reactionIndexs[e2->getId()]);
        }
        flowNum.push_back(bondIndexs[bnd]);
        flows.push_back("-f_" + bnd);
        efforts.push_back("e_" + bnd);
      }

      setDirac_[ie->getId()] = flowNum;
      SymbolicMatrix fi(flows.size(), 1);
      for (int f = 0; f < flows.size(); f++) {
        fi(f, 0) = SymEngine::parse(flows.at(f));
      }
      SymbolicMatrix ei(efforts.size(), 1);
      for (int e = 0; e < efforts.size(); e++) {
        ei(e, 0) = SymEngine::parse(efforts.at(e));
      }
      dirac[1] = fi;
      dirac[3] = ei;
      setDirac[ie->getId()] = dirac;
#ifdef DEBUG_PHS
      std::cout << elementIndexs[ie->getId()] << " " << ie->getName() << "\n";
      for (auto &c : flowNum)
        std::cout << "\t" << c << std::endl;
      std::cout << "\tfi\n" << fi << std::endl;
      std::cout << "\n\tei\n" << ei << std::endl;
#endif
    }

    // for all 0-elements
    for (auto &ie : setV0) {
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bin = giveBondsIncident(ie, setB);
      std::vector<std::string> fMat;
      std::vector<std::string> EMat;
      // Define F-matrix
      int n = bin.size();
      SymbolicMatrix Fi(n, n);
      for (int i = 0; i < n; i++) {
        Fi(0, i) = symOne;
      }
#ifdef DEBUG_PHS
      std::cout << " Zero junction " << n << std::endl;
      std::cout << "F matrix " << std::endl;
      std::cout << Fi << std::endl;
#endif
      // Define E-matrix
      SymbolicMatrix Ei(n, n);
      for (int i = 1; i < n; i++) {
        Ei(i, 0) = symOne;
        for (int j = 1; j < n; j++) {
          if (i == j)
            Ei(i, j) = symNegOne;
        }
      }
#ifdef DEBUG_PHS
      std::cout << "E matrix " << std::endl;
      std::cout << Ei << std::endl;
#endif
      setDirac[ie->getId()][0] = Fi;
      setDirac[ie->getId()][2] = Ei;
    }
    // for all 1-elements
    for (auto &ie : setV1) {
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bondsIn;
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bondsOut;
      std::tie(bondsIn, bondsOut) = giveBondsInOut(ie, setB, elementIndexs);
      auto mIni = bondsIn.size();   // Number of ingoing bonds in element ie
      auto mOuti = bondsOut.size(); // Number of outgoing bonds of element ie
      auto mi = mIni + mOuti;

      SymbolicMatrix matTi(mi, mi);
      for (int i = 0; i < mi; i++) {
        for (int j = 0; j < mi; j++) {
          if (i == j) {
            if (i < mIni) {
              matTi(i, j) = symOne;
            } else {
              matTi(i, j) = symNegOne;
            }
          }
        }
      }

      SymbolicMatrix thetaI(mi, mi);
      thetaI.fill(symZero);
      for (int i = 1; i < mi; i++) {
        thetaI(i, 0) = symOne;
        for (int j = 1; j < mi; j++) {
          if (i == j)
            thetaI(i, j) = symNegOne;
        }
      }
      SymbolicMatrix phiI(mi, mi);
      for (int i = 0; i < mi; i++) {
        phiI(0, i) = symOne;
      }
      // F Matrix
      auto Fi = thetaI * matTi;
      // E Matrix
      auto Ei = phiI * matTi;
#ifdef DEBUG_PHS
      std::cout << " One junction " << ie->getId() << "\t" << mi << " = "
                << mIni << " + " << mOuti << std::endl;
      std::cout << Fi << std::endl;
      std::cout << " Ei " << std::endl;
      std::cout << Ei << std::endl;
#endif
      setDirac[ie->getId()][0] = Fi;
      setDirac[ie->getId()][2] = Ei;
    }

    // for all Reaction-elements
    // Add in the 1 junction Dirac structure
    for (auto &ie : setVRe) {
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bondsIn;
      std::vector<std::tuple<RCPLIB::RCP<BGElement>, RCPLIB::RCP<BGElement>>>
          bondsOut;
      std::tie(bondsIn, bondsOut) = giveBondsInOut(ie, setB, elementIndexs);
      auto mIni = bondsIn.size();   // Number of ingoing bonds in element ie
      auto mOuti = bondsOut.size(); // Number of outgoing bonds of element ie
      auto mi = mIni + mOuti;

      SymbolicMatrix matTi(mi, mi);
      for (int i = 0; i < mi; i++) {
        for (int j = 0; j < mi; j++) {
          if (i == j) {
            if (i < mIni) {
              matTi(i, j) = symOne;
            } else {
              matTi(i, j) = symNegOne;
            }
          }
        }
      }

      SymbolicMatrix thetaI(mi, mi);
      thetaI.fill(symZero);
      for (int i = 1; i < mi; i++) {
        thetaI(i, 0) = symOne;
        for (int j = 1; j < mi; j++) {
          if (i == j)
            thetaI(i, j) = symNegOne;
        }
      }
      SymbolicMatrix phiI(mi, mi);
      for (int i = 0; i < mi; i++) {
        phiI(0, i) = symOne;
      }
      // F Matrix
      auto Fi = thetaI * matTi;
      // E Matrix
      auto Ei = phiI * matTi;
#ifdef DEBUG_PHS
      std::cout << " VRe 1 junction " << ie->getId() << "\t" << mi << " = "
                << mIni << " + " << mOuti << std::endl;
      std::cout << Fi << std::endl;
      std::cout << " Ei " << std::endl;
      std::cout << Ei << std::endl;
#endif
      setDirac[ie->getId()][0] = Fi;
      setDirac[ie->getId()][2] = Ei;
    }

    // for all Stoichiometry-elements
    // Constitutive equations are of the form f_0/r0 - f_1/r1 = 0; r0 * e_0 + r1
    // * e_1 = 0; F matrix is of the form {{1/r0,-1/r1}, {0,0}}, E matrix is of
    // the form {{0,0},{r0,r1}}.

    for (auto &ie : setVStoic) {
      int numStates = ie->getNumStates();
      auto values = ie->values();
      std::string pname1 = std::get<1>(values[numStates])->name + "_" +
                           std::to_string(elementIndexs[ie->getId()]);
      std::string pname2 = std::get<1>(values[numStates + 1])->name + "_" +
                           std::to_string(elementIndexs[ie->getId()]);

      SymbolicMatrix Fi(2, 2);
      Fi(0, 0) = SymEngine::parse(pname1);
      Fi(0, 1) = SymEngine::parse("-" + pname2);

      SymbolicMatrix Ei(2, 2);
      Ei(1, 0) = SymEngine::parse(pname1);
      Ei(1, 1) = SymEngine::parse(pname2);

      setDirac[ie->getId()][0] = Fi;
      setDirac[ie->getId()][2] = Ei;
    }

    // for all TF-elements
    for (auto &ie : setVTF) {
      int numStates = ie->getNumStates();
      auto values = ie->values();
      std::string pname = std::get<1>(values[numStates])->name + "_" +
                          std::to_string(elementIndexs[ie->getId()]);

      SymbolicMatrix Fi(2, 2);
      Fi(0, 0) = symOne;
      Fi(0, 1) = SymEngine::parse(pname);

      SymbolicMatrix Ei(2, 2);
      Ei(1, 1) = symOne;
      Ei(1, 0) = SymEngine::parse("-" + pname);

      setDirac[ie->getId()][0] = Fi;
      setDirac[ie->getId()][2] = Ei;
    }

    // for all GY-elements
    for (auto &ie : setVGY) {
      int numStates = ie->getNumStates();
      auto values = ie->values();
      std::string pname = std::get<1>(values[numStates])->name + "_" +
                          std::to_string(elementIndexs[ie->getId()]);

      SymbolicMatrix Fi(2, 2);

      Fi(0, 1) = SymEngine::parse(pname);
      Fi(1, 0) = SymEngine::parse("-" + pname);

      SymbolicMatrix Ei(2, 2);
      Ei(0, 0) = symOne;
      Ei(1, 1) = symOne;

      setDirac[ie->getId()][0] = Fi;
      setDirac[ie->getId()][2] = Ei;
    }

    // Permute the dirac matrices (rows) to get the bonds in the same order as
    // in setB (external elements followed by internals elements), using just
    // flows

    std::vector<int> setBEidx;
    std::vector<int> setBIidx;

    for (auto &b : setBE) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      setBEidx.push_back(bondIndexs[bnd]);
    }

    for (auto &b : setBI) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      setBIidx.push_back(bondIndexs[bnd]);
    }

#ifdef DEBUG_PHS
    // Inverse map of bond id to bond name, used for debugging
    std::map<int, std::string> indexBonds;
    for (auto &m : bondIndexs) {
      indexBonds[m.second] = m.first;
    }
#endif

    for (auto &d : setDirac_) {
      auto eid = d.first;
      std::vector<int> currentVectori = d.second;
#ifdef DEBUG_PHS
      std::cout << eid << " \nCurrent order :";
      for (auto &c : currentVectori) {
        std::cout << indexBonds[c] << "(" << c << ") ";
      }

      std::cout << std::endl;
      std::cout << " Exterior bonds ";
      for (auto &c : setBEidx) {
        std::cout << indexBonds[c] << "(" << c << ") ";
      }
      std::cout << "\n Interior bonds ";
      for (auto &c : setBIidx) {
        std::cout << indexBonds[c] << "(" << c << ") ";
      }
      std::cout << std::endl << "Wanted order :";
#endif
      std::vector<int> wantedVectori;
      // Maintain the order in the current vector, only moving interior elements
      // to the end
      for (auto &c : setBEidx) {
        auto ix = std::find(currentVectori.begin(), currentVectori.end(), c);
        if (ix != currentVectori.end())
          wantedVectori.push_back(c);
      }
      for (auto &c : setBIidx) {
        auto ix = std::find(currentVectori.begin(), currentVectori.end(), c);
        if (ix != currentVectori.end())
          wantedVectori.push_back(c);
      }

#ifdef DEBUG_PHS
      for (auto &cx : wantedVectori) {
        std::cout << indexBonds[cx] << "(" << cx << ") ";
      }
      std::cout << std::endl;
#endif
      auto matT = SymbolicMatrix(currentVectori.size(), currentVectori.size());
      matT.fill(symZero);
      for (int ci = 0; ci < currentVectori.size(); ci++) {
        auto cvi = currentVectori[ci];
        auto it = std::find(wantedVectori.begin(), wantedVectori.end(), cvi);
        // If element was found
        if (it != wantedVectori.end()) {
          int mi = it - wantedVectori.begin();
          matT(mi, ci) = symOne;
        }
      }

#ifdef DEBUG_PHS
      std::cout << "Transformation matrix \n" << matT << std::endl;
#endif

      auto Fi = setDirac[eid][0];
      auto fi = setDirac[eid][1];
      auto Ei = setDirac[eid][2];
      auto ei = setDirac[eid][3];

#ifdef DEBUG_PHS
      std::cout << "Fi \n" << Fi << std::endl;
      std::cout << "fi \n" << fi << std::endl;
      std::cout << "Ei \n" << Fi << std::endl;
      std::cout << "ei \n" << Ei << std::endl;
#endif

      setDirac[eid][0] = Fi * matT.transpose();
      setDirac[eid][1] = (matT * fi);

      setDirac[eid][2] = Ei * matT.transpose();
      setDirac[eid][3] = matT * ei;

#ifdef DEBUG_PHS
      std::cout << "Mat T " << compIDMap[eid]->getName() << std::endl;
      std::cout << matT << "\n" << std::endl;

      std::cout << "Flow  \n"
                << Fi << "\n After \n"
                << setDirac[eid][0] << std::endl;

      std::cout << "Flow vec \n"
                << fi << "\n After \n"
                << setDirac[eid][1] << std::endl;

      std::cout << "Effort  \n"
                << Ei << "\n After \n"
                << setDirac[eid][2] << std::endl;

      std::cout << "Effort vec \n"
                << ei << "\n After \n"
                << setDirac[eid][3] << std::endl;
#endif
    }

    // Compute D_IC
    int sBIsize = setBI.size();
    SymbolicMatrix vecFK(sBIsize * 2, 1);
    SymbolicMatrix vecEK(sBIsize * 2, 1);
    int vix = 0;
    for (auto &b : setBI) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance

        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      vecFK(vix, 0) = SymEngine::parse("f_" + bnd);
      vecFK(vix + sBIsize, 0) = SymEngine::parse("-f_" + bnd);
      vecEK(vix, 0) = SymEngine::parse("e_" + bnd);
      vecEK(vix + sBIsize, 0) = SymEngine::parse("e_" + bnd);

      vix++;
    }

    std::vector<std::string> setDiracKeys;
    // Ensure the order 0, 1, Re, Stoic, TF, GY
    for (auto &ie : setV0) {
      setDiracKeys.push_back(ie->getId());
    }

    for (auto &ie : setV1) {
      setDiracKeys.push_back(ie->getId());
    }

    for (auto &ie : setVRe) {
      setDiracKeys.push_back(ie->getId());
    }

    for (auto &ie : setVStoic) {
      setDiracKeys.push_back(ie->getId());
    }

    for (auto &ie : setVTF) {
      setDiracKeys.push_back(ie->getId());
    }

    for (auto &ie : setVGY) {
      setDiracKeys.push_back(ie->getId());
    }

    SymbolicMatrix vecFIC(setDirac[setDiracKeys.at(0)][1]);

    for (int i = 1; i < setDiracKeys.size(); i++) {
      SymbolicMatrix v = setDirac[setDiracKeys.at(i)][1];
      SymbolicMatrix cv(vecFIC.rows() + v.rows(), 1);
      cv << vecFIC, v;
      vecFIC = cv;
    }

    // Get the permutation matrix for rearranging rows to follow the order
    // in which the flows appear in setBI

#ifdef DEBUG_PHS
    std::cout << vecFK << std::endl << std::endl;
    std::cout << vecFIC << std::endl << std::endl;
#endif

    std::vector<int> forder;
    std::vector<int> sorder;
    for (int ci = 0; ci < sBIsize * 2; ci++) {
      int it = getEigenVectorIndex(vecFK(ci, 0), vecFIC);
      if (it != -1) {
        forder.push_back(it);
        sorder.push_back(it);
      } else {
        // Error report
        logDebug("Error finding the corresponding eigen vector !!");
      }
    }
    // Get the sorted index
    std::sort(sorder.begin(), sorder.end());
    std::vector<int> sindex(forder.size());
    for (int ci = 0; ci < sBIsize * 2; ci++) {
      auto v = forder[ci];
      auto it = std::find(sorder.begin(), sorder.end(), v) - sorder.begin();
      sindex[ci] = it;
    }
    auto vmatT = SymbolicMatrix(sBIsize * 2, sBIsize * 2);
    vmatT.fill(symZero);

    for (int ci = 0; ci < sBIsize * 2; ci++) {
      vmatT(sindex[ci], ci) = symOne;
    }

    SymbolicMatrix matFIC(2 * numMI, 2 * numMI);
    SymbolicMatrix matEIC(2 * numMI, 2 * numMI);
    matFIC.fill(symZero);
    matEIC.fill(symZero);
    auto &I1 = matFIC.block(0, 0, numMI, numMI);
    auto &I2 = matFIC.block(0, numMI, numMI, numMI);
    auto &E1 = matEIC.block(numMI, 0, numMI, numMI);
    auto &E2 = matEIC.block(numMI, numMI, numMI, numMI);

    for (int i = 0; i < numMI; i++) {
      I1(i, i) = symOne;
      I2(i, i) = symOne;
      E1(i, i) = symOne;
      E2(i, i) = symNegOne;
    }

    auto DIC_Em = matEIC * vmatT.transpose();
    auto DIC_Fm = matFIC * vmatT.transpose();
    auto DIC_Ev = vmatT * vecEK;
    auto DIC_Fv = vmatT * vecFK;

    std::map<std::string, int> setNumME;
    std::map<std::string, int> setNumMI;
    std::map<std::string, int> setNumM;
    for (auto &el : setDiracKeys) {
      auto elem = compIDMap[el];
      setNumME[el] = giveBondsIncident(elem, setBE).size();
      setNumMI[el] = giveBondsIncident(elem, setBI).size();
      setNumM[el] = setNumMI[el] + setNumME[el];
#ifdef DEBUG_PHS
      std::cout << elem->getName() << " E " << setNumME[el] << " I "
                << setNumMI[el] << std::endl;
#endif
    }

    // Initialize matrices F_IC(i) and E_IC(i) of interconnection Dirac
    // structure

    std::map<std::string, SymbolicMatrix> matFICi;
    std::map<std::string, SymbolicMatrix> matEICi;
    unsigned int im = 0;

    for (auto &c : setDiracKeys) {
#ifdef DEBUG_PHS
      std::cout << "Col indexes " << im << "\t" << im + setNumMI[c]
                << std::endl;
#endif
      std::vector<unsigned int> ind = {im};
      if (setNumMI[c] > 1) {
        for (int ci = im + 1; ci < im + setNumMI[c]; ci++)
          ind.push_back(ci);
      }
      matFICi[c] = DIC_Fm(Eigen::placeholders::all, ind);
      matEICi[c] = DIC_Em(Eigen::placeholders::all, ind);
#ifdef DEBUG_PHS
      std::cout << compIDMap[c]->getName() << "\nFIC\n"
                << matFICi[c] << "\nEIC\n"
                << matEICi[c] << std::endl;
#endif
      im += setNumMI[c];
    }

    // matrices M^T (i)
    std::map<std::string, SymbolicMatrix> matMTi;
    int numCols = 0; // The number of columns for Gamma matrix
    for (auto &c : setDiracKeys) {
      auto sdf = setDirac[c][0].transpose();
      auto sde = setDirac[c][2].transpose();
      std::vector<unsigned int> ind;
      for (int ci = 0; ci < setNumME[c]; ci++) {
        ind.push_back(ci);
      }

      setDirac[c][4] = sdf(ind, Eigen::placeholders::all);
      setDirac[c][6] = sde(ind, Eigen::placeholders::all);
      ind.clear();
      for (int ci = setNumME[c]; ci < sdf.rows(); ci++) {
        ind.push_back(ci);
      }
      setDirac[c][5] = sdf(ind, Eigen::placeholders::all);
      setDirac[c][7] = sde(ind, Eigen::placeholders::all);
      matMTi[c] = matFICi[c] * setDirac[c][7] + matEICi[c] * setDirac[c][5];

#ifdef DEBUG_PHS
      std::cout << compIDMap[c]->getName() << "F \n" << sdf << std::endl;
      std::cout << compIDMap[c]->getName() << " F first \n"
                << setDirac[c][4] << std::endl;
      std::cout << compIDMap[c]->getName() << " F reminder\n"
                << setDirac[c][5] << std::endl;
      std::cout << compIDMap[c]->getName() << "MatMTi\n"
                << matMTi[c] << std::endl;
#endif

      numCols += matMTi[c].cols();
    }
    SymbolicMatrix matMT(2 * numMI, numCols);
    matMT.fill(symZero);
    int colCounter = 0;
    for (auto &c : setDiracKeys) {
      int nCols = matMTi[c].cols();
      matMT.block(0, colCounter, matMT.rows(), nCols) = matMTi[c];
      colCounter += nCols;
    }
#ifdef DEBUG_PHS
    std::cout << "Gamma matrix \n" << matMT << std::endl;
#endif
    // Find the nullspace - following the sympy implementation
    SymEngine::DenseMatrix s_MT(matMT.rows(), matMT.cols());
    SymEngine::DenseMatrix reduced(matMT.rows(), matMT.cols());
    for (int r = 0; r < matMT.rows(); r++) {
      for (int c = 0; c < matMT.cols(); c++) {
        s_MT.set(r, c, matMT(r, c));
      }
    }

    SymEngine::vec_uint pivots;
    SymEngine::reduced_row_echelon_form(s_MT, reduced, pivots);

    SymEngine::vec_uint free_vars;
    for (int i = 0; i < matMT.cols(); i++) {
      auto it = std::find(pivots.begin(), pivots.end(), i);
      if (it == pivots.end()) {
        free_vars.push_back(i);
      }
    }
    std::vector<SymEngine::vec_basic> basis;
    for (auto &fv : free_vars) {
      SymEngine::vec_basic vec(matMT.cols(), symZero);
      vec[fv] = symOne;

      for (int piv_row = 0; piv_row < pivots.size(); piv_row++) {
        // Handle multiple overload of - operator, by casting 'vec' to
        // Expression
        SymEngine::Expression v1 = vec[pivots[piv_row]];
        SymEngine::Expression v2 = reduced.get(piv_row, fv);

        vec[pivots[piv_row]] = v1 - v2;
      }
      basis.push_back(vec);
    }
    // The order of free_variables is revered in the mathematica code/ paper
    // example So following it through
    int numBasis = basis.size();
    SymbolicMatrix matL(numBasis, matMT.cols());
    matL.fill(symZero);
    for (int i = 0; i < numBasis; i++) {
      SymEngine::vec_basic &vec = basis.back();
      for (int j = 0; j < vec.size(); j++) {
        matL(i, j) = vec[j];
      }
      basis.pop_back();
    }

    // Extract block matrices L (i) from L
    std::map<std::string, SymbolicMatrix> matLi;
    colCounter = 0;
    for (auto &c : setDiracKeys) {
      std::vector<unsigned int> ind;
      for (int ci = 0; ci < setNumM[c]; ci++) {
        ind.push_back(ci + colCounter);
      }
      matLi[c] = matL(Eigen::placeholders::all, ind);
      colCounter += setNumM[c];
    }

    std::map<int, SymbolicMatrix> dirac;
    for (auto &c : setDiracKeys) {
      if (setNumME[c] > 0) {

        auto dr1 = matLi[c] * (setDirac[c][4].transpose());
        auto dr3 = matLi[c] * (setDirac[c][6].transpose());
        std::vector<unsigned int> ind;
        for (int ci = 0; ci < setNumME[c]; ci++) {
          ind.push_back(ci);
        }
        auto dr2 = setDirac[c][1](ind, Eigen::placeholders::all);
        auto dr4 = setDirac[c][3](ind, Eigen::placeholders::all);

        if (dirac.find(1) != dirac.end()) {
          SymbolicMatrix newm(dirac[1].rows(), dirac[1].cols() + dr1.cols());
          newm << dirac[1], dr1;
          dirac[1] = newm;
        } else {
          dirac[1] = dr1;
        }
        if (dirac.find(3) != dirac.end()) {
          SymbolicMatrix newm(dirac[3].rows(), dirac[3].cols() + dr1.cols());
          newm << dirac[3], dr3;
          dirac[3] = newm;
        } else {
          dirac[3] = dr3;
        }
        if (dirac.find(2) != dirac.end()) {
          SymbolicMatrix newm(dirac[2].rows() + dr2.rows(), dirac[2].cols());
          newm << dirac[2], dr2;
          dirac[2] = newm;
        } else {
          dirac[2] = dr2;
        }
        if (dirac.find(4) != dirac.end()) {
          SymbolicMatrix newm(dirac[4].rows() + dr4.rows(), dirac[4].cols());

          newm << dirac[4], dr4;
          dirac[4] = newm;
        } else {
          dirac[4] = dr4;
        }
#ifdef DEBUG_PHS
        std::cout << compIDMap[c]->getName() << "\n"
                  << dirac[1] << "\n\n"
                  << dirac[2] << "\n\n"
                  << dirac[3] << "\n\n"
                  << dirac[4] << std::endl;
#endif
      }
    }

    // Compute permuation matrix T
    std::vector<std::string> dvec;
    SymbolicMatrix &no = dirac[2];
    for (int i = 0; i < no.rows(); i++) {
      std::string c = no(i, 0).get_basic()->__str__();
      auto loc = c.find('_');
      dvec.push_back(c.substr(loc + 1, c.size()));
    }

    SymbolicMatrix matT(dvec.size(), dvec.size());
    std::vector<std::string> wantedVector;
    std::vector<std::string> setBEbnd;
    std::vector<std::string> setBIbnd;
    std::map<std::string, std::string> bondElementMap;

    for (auto &b : setBE) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance

        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      bondElementMap[bnd] = e1->getName() + "->" + e2->getName();
      setBEbnd.push_back(bnd);
    }

    for (auto &b : setBI) {
      RCPLIB::RCP<BGElement> &e1 = std::get<0>(b);
      RCPLIB::RCP<BGElement> &e2 = std::get<1>(b);
      std::string bnd;
      if (e1->getId() != e2->getId()) { // Handle VRe elements which have both a
                                        // 1 junction and a resistance

        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(elementIndexs[e2->getId()]);
      } else {
        bnd = std::to_string(elementIndexs[e1->getId()]) + "_" +
              std::to_string(reactionIndexs[e2->getId()]);
      }
      bondElementMap[bnd] = e1->getName() + "->" + e2->getName();
      setBIbnd.push_back(bnd);
    }
    // Store Bond element mapping
    result["BondIdMap"] = bondElementMap;

    // Maintain the order in the current vector, only moving interior elements
    // to the end
    for (auto &c : setBEbnd) {
      auto ix = std::find(dvec.begin(), dvec.end(), c);
      if (ix != dvec.end())
        wantedVector.push_back(*ix);
    }
    for (auto &c : setBIbnd) {
      auto ix = std::find(dvec.begin(), dvec.end(), c);
      if (ix != dvec.end())
        wantedVector.push_back(*ix);
    }

    auto matTd = SymbolicMatrix(dvec.size(), dvec.size());
    matTd.fill(symZero);
    for (int ci = 0; ci < dvec.size(); ci++) {
      auto cvi = dvec[ci];
      auto it = std::find(wantedVector.begin(), wantedVector.end(), cvi);
      // If element was found
      if (it != wantedVector.end()) {
        int mi = it - wantedVector.begin();
        matTd(mi, ci) = symOne;
      }
    }

    auto dr1 = dirac[1] * matTd.transpose();
    auto dr3 = dirac[3] * matTd.transpose();
    auto dr2 = matTd * dirac[2];
    auto dr4 = matTd * dirac[4];

    dirac[1] = dr1;
    dirac[2] = dr2;
    dirac[3] = dr3;
    dirac[4] = dr4;
#ifdef DEBUG_PHS
    std::cout << "Dirac 1\n"
              << dr1 << std::endl
              << "Dirac 2\n"
              << dr2 << std::endl
              << "Dirac 3\n"
              << dr3 << std::endl
              << "Dirac 4\n"
              << dr4 << std::endl;
#endif
    // Store implicit DS
    nlohmann::json ids;
    ids["F"] = to_json(dirac[1]);
    ids["f"] = to_json(dirac[2]);
    ids["E"] = to_json(dirac[3]);
    ids["e"] = to_json(dirac[4]);
    result["ImplicitDS"] = ids;

#ifdef CHECK_DEPENDEND_SOURCES
    // Expensive tests - ignore if possible
    // Test Assumption rank(E_Sf | F_Se) = N(nSf + nSe)
    if (numNS > 0) {
      std::vector<int> ind;
      int offset = dirac[3].cols() - numNS;
      for (int i = 0; i < numNS - numNSe; i++) {
        ind.push_back(offset + i);
      }
      auto E_Sf = dirac[3](Eigen::placeholders::all, ind);
      ind.clear();
      offset = dirac[1].cols() - numNSe;
      for (int i = 0; i < numNSe; i++) {
        ind.push_back(offset + i);
      }
      auto F_Se = SymbolicMatrix(0, 0);
      Eigen::MatrixX<float> dm;

      int rank = 0;
      if (ind.size() > 0) {
        F_Se = dirac[1](Eigen::placeholders::all, ind);
        SymbolicMatrix ss(E_Sf.rows(), E_Sf.cols() + F_Se.cols());
        ss << E_Sf, F_Se;
        dm = Eigen::MatrixX<float>(ss.rows(), ss.cols());
        for (int r = 0; r < ss.rows(); r++) {
          for (int c = 0; c < ss.cols(); c++) {
            if (!(ss(r, c) == symZero)) {
              dm(r, c) = 1;
            }
          }
        }
        Eigen::FullPivLU<Eigen::MatrixX<float>> lu_decomp(dm);
        rank = lu_decomp.rank();

      } else {
        dm = Eigen::MatrixX<float>(E_Sf.rows(), E_Sf.cols());
        for (int r = 0; r < E_Sf.rows(); r++) {
          for (int c = 0; c < E_Sf.cols(); c++) {
            if (!(E_Sf(r, c) == symZero)) {
              dm(r, c) = 1;
            }
          }
        }
        Eigen::FullPivLU<Eigen::MatrixX<float>> lu_decomp(dm);
        rank = lu_decomp.rank();
      }

      if (rank != numNS)
        warnings.push_back(
            "The assumption rank(E_Sf | F_Se) = N(nSf + nSe) is not fulfiled. "
            "The bond Graph contains dependent sources.");
      // throw BGException(
      //     "The assumption rank(E_Sf | F_Se) = N(nSf + nSe) is not
      //     fulfiled. " "The bond Graph contains dependent sources.");

      // Test Assumption rank(F_C | E_Sf | F_Se) = N(nC + nSf + nSe)

      ind.clear();
      for (int i = 0; i < numNC; i++) {
        ind.push_back(i);
      }
      if (ind.size() > 0) {
        auto F_C = dirac[1](Eigen::placeholders::all, ind);
        Eigen::MatrixX<float> fcm(F_C.rows(), F_C.cols() + dm.cols());
        for (int r = 0; r < F_C.rows(); r++) {
          for (int c = 0; c < F_C.cols(); c++) {
            if (!(F_C(r, c) == symZero)) {
              fcm(r, c) = 1;
            }
          }
        }
        fcm.block(0, F_C.cols(), dm.rows(), dm.cols()) = dm;
        Eigen::FullPivLU<Eigen::MatrixX<float>> lu_decomp(fcm);
        rank = lu_decomp.rank();
      }

      if (rank != (numNC + numNS))
        // throw BGException(
        //     "The assumption rank(F_C | E_Sf | F_Se) = N(nC + nSf + nSe) is
        //     not \
        //   fulfiled. The bond graph contains dependent storages or storages \
        //   determined by sources.");
        warnings.push_back(
            "The assumption rank(F_C | E_Sf | F_Se) = N(nC + nSf + nSe) is not "
            "fulfiled. The bond Graph contains dependent sources.");
    }
#endif

    // Determine splitting of resistive elements to R1 and R2
    std::vector<unsigned int> numNEAcc;
    std::vector<int> posR1;
    std::vector<int> posR2;
    std::vector<int> ind;
    if (numNR > 0) {
      // Calculate accumulated number of elements {nC, nC + nR, ...}
      numNEAcc.push_back(numNC);
      numNEAcc.push_back(numNC + numNR);
      numNEAcc.push_back(numNC + numNR + numNSf);
      numNEAcc.push_back(numNC + numNR + numNSf + numNSe);

      // Concatenate block matrices to matPartZ1 = {F_C, E_Sf, F_Se}
      ind.clear();
      for (int i = 0; i < numNEAcc[0]; i++)
        ind.push_back(i);
      auto F_C = dirac[1](Eigen::placeholders::all, ind);
      ind.clear();
      for (int i = numNEAcc[1]; i < numNEAcc[2]; i++)
        ind.push_back(i);
      auto E_Sf = dirac[3](Eigen::placeholders::all, ind);
      ind.clear();
      for (int i = numNEAcc[2]; i < numNEAcc[3]; i++)
        ind.push_back(i);
      auto F_Se = dirac[1](Eigen::placeholders::all, ind);

      SymbolicMatrix matPartZ1(F_C.rows(),
                               F_C.cols() + E_Sf.cols() + F_Se.cols());
      matPartZ1 << F_C, E_Sf, F_Se;
#ifdef DEBUG_PHS
      std::cout << " matPartZ1 " << std::endl;
      std::cout << matPartZ1 << std::endl;
#endif
      int rank = getMatrixRank(matPartZ1);

      // Determine splitting of elements R to R1 and R2
      for (int iR = numNEAcc[0]; iR < numNEAcc[1]; iR++) {
        ind.clear();
        ind.push_back(iR);
        auto matPart2 = dirac[1](Eigen::placeholders::all, ind);

        SymbolicMatrix augPart(matPartZ1.rows(),
                               matPartZ1.cols() + matPart2.cols());
        augPart << matPartZ1, matPart2;
        auto arank = getMatrixRank(augPart);
#ifdef DEBUG_PHS
        std::cout << iR << " " << augPart << std::endl;
        std::cout << "Rank " << rank << " " << arank << std::endl;
#endif
        if (arank == rank + 1) {
          matPartZ1 = augPart;
          rank++;
          posR1.push_back(iR);
        } else {
          posR2.push_back(iR);
        }
      }

      // Compute permuation matrix T that is an identity matrix except for
      // entries belonging to resistive elements
      std::vector<std::string> vecFAbs; // Absolute value of vector of f
      for (int i = 0; i < dirac[2].rows(); i++) {
        std::string v = dirac[2](i, 0).get_basic()->__str__();
        auto loc = v.find('_');
        vecFAbs.push_back(v.substr(loc + 1, v.size()));
      }
      std::vector<std::string> currentVector;
      for (int i = numNEAcc[0]; i < numNEAcc[1]; i++) {
        currentVector.push_back(vecFAbs[i]);
      }
      std::vector<std::string> wantedVector;
      for (auto &c : posR1) {
        wantedVector.push_back(vecFAbs[c]);
      }
      for (auto &c : posR2) {
        wantedVector.push_back(vecFAbs[c]);
      }

      auto matTR = SymbolicMatrix(currentVector.size(), currentVector.size());
      matTR.fill(symZero);
      for (int ci = 0; ci < currentVector.size(); ci++) {
        auto cvi = currentVector[ci];
        auto it = std::find(wantedVector.begin(), wantedVector.end(), cvi);
        // If element was found
        if (it != wantedVector.end()) {
          int mi = it - wantedVector.begin();
          matTR(mi, ci) = symOne;
        }
      }

      SymbolicMatrix matT(numNEAcc[3], numNEAcc[3]);
      matT.fill(symZero);
      for (int i = 0; i < numNEAcc[3]; i++)
        matT(i, i) = symOne;

      matT.block(numNEAcc[0], numNEAcc[0], matTR.rows(), matTR.cols()) = matTR;

      auto dr1 = dirac[1] * matT.transpose();
      auto dr3 = dirac[3] * matT.transpose();
      auto dr2 = matT * dirac[2];
      auto dr4 = matT * dirac[4];
      dirac[1] = dr1;
      dirac[2] = dr2;
      dirac[3] = dr3;
      dirac[4] = dr4;
    }
#ifdef DEBUG_PHS
    std::cout << "Implicit dirac structure F f E e" << std::endl;
    std::cout << dirac[1] << std::endl << dirac[2] << std::endl;
    std::cout << dirac[3] << std::endl << dirac[4] << std::endl;
#endif
    std::map<int, SymbolicMatrix> diracEx;
    numNEAcc.clear();
    // Calculate accumulated number of elements {nC, nC + nR1, ...}
    numNEAcc.push_back(setVC.size() + setVConc.size());
    numNEAcc.push_back(setVC.size() + setVConc.size() + posR1.size());
    numNEAcc.push_back(setVC.size() + setVConc.size() + posR1.size() +
                       posR2.size());
    numNEAcc.push_back(setVC.size() + setVConc.size() + posR1.size() +
                       posR2.size() + setVSf.size());
    numNEAcc.push_back(setVC.size() + setVConc.size() + posR1.size() +
                       posR2.size() + setVSf.size() + setVSe.size());

    // Concatenate vector veyY = {f_C, f_R1, e_R2, e_Sf, f_Se}
    ind.clear();
    for (int i = 0; i < numNEAcc[0]; i++)
      ind.push_back(i);
    auto f_C = dirac[2](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[0]; i < numNEAcc[1]; i++)
      ind.push_back(i);
    auto f_R1 = dirac[2](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[1]; i < numNEAcc[2]; i++)
      ind.push_back(i);
    auto e_R2 = dirac[4](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[2]; i < numNEAcc[3]; i++)
      ind.push_back(i);
    auto e_Sf = dirac[4](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[3]; i < numNEAcc[4]; i++)
      ind.push_back(i);
    auto f_Se = dirac[2](ind, Eigen::placeholders::all);

    SymbolicMatrix dex1(
        f_C.rows() + f_R1.rows() + e_R2.rows() + e_Sf.rows() + f_Se.rows(), 1);
    dex1 << f_C, f_R1, e_R2, e_Sf, f_Se;
    diracEx[1] = dex1;

    ind.clear();
    for (int i = 0; i < numNEAcc[0]; i++)
      ind.push_back(i);
    auto e_C = dirac[4](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[0]; i < numNEAcc[1]; i++)
      ind.push_back(i);
    auto e_R1 = dirac[4](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[1]; i < numNEAcc[2]; i++)
      ind.push_back(i);
    auto f_R2 = dirac[2](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[2]; i < numNEAcc[3]; i++)
      ind.push_back(i);
    auto f_Sf = dirac[2](ind, Eigen::placeholders::all);
    ind.clear();
    for (int i = numNEAcc[3]; i < numNEAcc[4]; i++)
      ind.push_back(i);
    auto e_Se = dirac[4](ind, Eigen::placeholders::all);

    SymbolicMatrix dex3(
        e_C.rows() + e_R1.rows() + f_R2.rows() + f_Sf.rows() + e_Se.rows(), 1);
    dex3 << e_C, e_R1, f_R2, f_Sf, e_Se;
    diracEx[3] = dex3;

    // Concatenate block matrices to matZ1 = {F_C, F_R1, E_R2, E_Sf, F_Se}
#ifdef DEBUG_PHS
    std::cout << "posR1 " << posR1.size() << std::endl;
    std::cout << "posR2 " << posR2.size() << std::endl;
    for (auto &c : numNEAcc) {
      std::cout << c << " ";
    }
    std::cout << std::endl;
#endif

    ind.clear();
    for (int i = 0; i < numNEAcc[0]; i++)
      ind.push_back(i);
    auto F_C = dirac[1](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[0]; i < numNEAcc[1]; i++)
      ind.push_back(i);
    auto F_R1 = dirac[1](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[1]; i < numNEAcc[2]; i++)
      ind.push_back(i);
    auto E_R2 = dirac[3](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[2]; i < numNEAcc[3]; i++)
      ind.push_back(i);
    auto E_Sf = dirac[3](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[3]; i < numNEAcc[4]; i++)
      ind.push_back(i);
    auto F_Se = dirac[1](Eigen::placeholders::all, ind);
#ifdef DEBUG_PHS
    std::cout << "FC\n" << F_C << std::endl;
    std::cout << "FR1\n" << F_R1 << std::endl;
    std::cout << "ER2\n" << E_R2 << std::endl;
    std::cout << "ESf\n" << E_Sf << std::endl;
    std::cout << "FSe\n" << F_Se << std::endl;
#endif
    SymbolicMatrix matZ1(F_C.rows(), F_C.cols() + F_R1.cols() + E_R2.cols() +
                                         E_Sf.cols() + F_Se.cols());
    matZ1 << F_C, F_R1, E_R2, E_Sf, F_Se;
    // Concatenate block matrices to matZ2 = {E_C, E_R1, F_R2, F_Sf, E_Se}
    ind.clear();
    for (int i = 0; i < numNEAcc[0]; i++)
      ind.push_back(i);
    auto E_C = dirac[3](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[0]; i < numNEAcc[1]; i++)
      ind.push_back(i);
    auto E_R1 = dirac[3](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[1]; i < numNEAcc[2]; i++)
      ind.push_back(i);
    auto F_R2 = dirac[1](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[2]; i < numNEAcc[3]; i++)
      ind.push_back(i);
    auto F_Sf = dirac[1](Eigen::placeholders::all, ind);
    ind.clear();
    for (int i = numNEAcc[3]; i < numNEAcc[4]; i++)
      ind.push_back(i);
    auto E_Se = dirac[3](Eigen::placeholders::all, ind);

    SymbolicMatrix matZ2(E_C.rows(), E_C.cols() + E_R1.cols() + F_R2.cols() +
                                         F_Sf.cols() + E_Se.cols());
#ifdef DEBUG_PHS
    std::cout << "EC\n" << E_C << std::endl;
    std::cout << "ER1\n" << E_R1 << std::endl;
    std::cout << "FR2\n" << F_R2 << std::endl;
    std::cout << "FSf\n" << F_Sf << std::endl;
    std::cout << "ESe\n" << E_Se << std::endl;
#endif
    matZ2 << E_C, E_R1, F_R2, F_Sf, E_Se;
#ifdef DEBUG_PHS
    std::cout << "\n" << matZ1 << std::endl << std::endl;
    std::cout << matZ2 << std::endl;
#endif
    SymEngine::DenseMatrix Z1(matZ1.rows(), matZ1.cols());
    SymEngine::DenseMatrix Z2(matZ2.rows(), matZ2.cols());
    for (int i = 0; i < matZ1.rows(); i++) {
      for (int j = 0; j < matZ1.cols(); j++) {
        Z1.set(i, j, matZ1(i, j));
      }
    }
    for (int i = 0; i < matZ2.rows(); i++) {
      for (int j = 0; j < matZ2.cols(); j++) {
        Z2.set(i, j, matZ2(i, j));
      }
    }

    SymEngine::DenseMatrix Z1inv(matZ1.rows(), matZ1.cols());
    SymEngine::DenseMatrix Z1invProdZ2(matZ1.rows(), matZ1.cols());
    try {
      Z1.inv(Z1inv);
      Z1inv.mul_matrix(Z2, Z1invProdZ2);
    } catch (SymEngine::SymEngineException
                 &e) { // Handle rank deficient Z1 matrix, when inv will fail
      Z1.transpose(Z1inv);
      Z1inv.mul_matrix(Z2, Z1invProdZ2);
    }
    SymbolicMatrix drex3(matZ1.rows(), matZ1.cols());
    for (int i = 0; i < matZ1.rows(); i++) {
      for (int j = 0; j < matZ1.cols(); j++) {
        SymEngine::Expression expr = Z1invProdZ2.get(i, j);
        drex3(i, j) = SymEngine::simplify(-expr);
      }
    }
    diracEx[2] = drex3;

    // Calculate number of bond graph elements
    auto numC = numNEAcc[0];               // Number of storages
    auto numR = numNEAcc[2] - numNEAcc[0]; // Number of resistors
    auto numS = numNEAcc[4] - numNEAcc[2]; // Number of sources

    result["number of storages"] = numC;
    result["number of resistances"] = numR;
    result["number of sources"] = numS;

    SymbolicMatrix zCC = drex3.block(0, 0, numC, numC);
    SymbolicMatrix zCR = -drex3.block(0, numC, numC, numR);
    SymbolicMatrix zCP = -drex3.block(0, drex3.cols() - numS, numC, numS);
    SymbolicMatrix zRR = drex3.block(numC, numC, numR, numR);
    SymbolicMatrix zRP = -drex3.block(numC, drex3.cols() - numS, numR, numS);
    SymbolicMatrix zPP =
        drex3.block(drex3.rows() - numS, drex3.cols() - numS, numS, numS);

#ifdef DEBUG_PHS
    std::cout << " zCC \n" << zCC << std::endl;
    std::cout << " zCR \n" << zCR << std::endl;
    std::cout << " zCP \n" << zCP << std::endl;
    std::cout << " zRR \n" << zRR << std::endl;
    std::cout << " zRP \n" << zRP << std::endl;
    std::cout << " zPP \n" << zPP << std::endl;
#endif

    SymbolicMatrix matA;
    SymbolicMatrix matK;
    SymbolicMatrix matB;
    SymbolicMatrix matJ;
    SymbolicMatrix matR;
    SymbolicMatrix matG;
    SymbolicMatrix matP;
    SymbolicMatrix matPT;
    SymbolicMatrix matM;
    SymbolicMatrix matS;
    SymbolicMatrix vecY;
    SymbolicMatrix vecU;

    if (numR > 0) {
      // TODO Test for multiple R
      // Get the resistances
      SymEngine::vec_basic resistance;
      Eigen::MatrixX<number> matD2;

      if (numNR > 1) {
        matD2 = SymbolicMatrix(numNR, 2);
        matD2.fill(symZero);
        int ix = 0;

        for (auto &ie : setVR) {
          int numStates = ie->getNumStates();
          auto values = ie->values();
          std::string pname = std::get<1>(values[numStates])->name + "_" +
                              std::to_string(elementIndexs[ie->getId()]);
          resistance.push_back(SymEngine::parse(pname));
          matD2(ix++, 1) = SymEngine::parse(pname);
        }
        for (auto &ie : setVRe) {
          int numStates = ie->getNumStates();
          auto values = ie->values();
          std::string pname = std::get<1>(values[numStates])->name + "_" +
                              std::to_string(reactionIndexs[ie->getId()]);
          resistance.push_back(SymEngine::parse(pname));
          matD2(ix++, 1) = SymEngine::parse(pname);
        }

        matD2(0, 0) = matD2(0, 1);
        matD2(0, 1) = symZero;
      } else {
        matD2 = SymbolicMatrix(numNR, 1);
        for (auto &ie : setVR) {
          int numStates = ie->getNumStates();
          auto values = ie->values();
          std::string pname = std::get<1>(values[numStates])->name + "_" +
                              std::to_string(elementIndexs[ie->getId()]);
          resistance.push_back(SymEngine::parse(pname));
          matD2(0, 0) = SymEngine::parse(pname);
        }
        for (auto &ie : setVRe) {
          int numStates = ie->getNumStates();
          auto values = ie->values();
          std::string pname = std::get<1>(values[numStates])->name + "_" +
                              std::to_string(reactionIndexs[ie->getId()]);
          resistance.push_back(SymEngine::parse(pname));
          matD2(0, 0) = SymEngine::parse(pname);
        }
      }

      std::vector<std::string> vecfR_;
      std::vector<std::string> vecuR_;
      ind.clear();
      for (int i = numNC; i < numNC + numNR; i++) {
        ind.push_back(i);
        vecfR_.push_back("f_" + setBEbnd[i]);
      }
      auto drur = diracEx[3](ind, Eigen::placeholders::all);
      for (int i = 0; i < drur.rows(); i++) {
        std::string expr = drur(i, 0).get_basic()->__str__();
        if (expr.at(0) == '-') {
          expr = expr.substr(1, expr.size());
        }
        vecuR_.push_back(expr);
      }
      SymbolicMatrix vecfR(vecfR_.size(), 1);
      for (int i = 0; i < vecfR_.size(); i++) {
        vecfR(i, 0) = SymEngine::parse(vecfR_[i]);
      }
      SymbolicMatrix vecuR(vecuR_.size(), 1);
      for (int i = 0; i < vecuR_.size(); i++) {
        vecuR(i, 0) = SymEngine::parse(vecuR_[i]);
      }

      SymbolicMatrix matT(vecfR.size(), vecfR.size());
      matT.fill(symZero);
      for (int ci = 0; ci < vecfR.size(); ci++) {
        auto cvi = vecfR_[ci];
        auto it = std::find(vecuR_.begin(), vecuR_.end(), cvi);
        // If element was found
        if (it != vecuR_.end()) {
          int mi = it - vecuR_.begin();
          matT(mi, ci) = symOne;
        }
      }

      auto vecfR1R2 = matT * vecfR;
      SymbolicMatrix matD2R1R2 = matD2 * matT;

      SymEngine::DenseMatrix D2R1R2(matD2R1R2.rows(), matD2R1R2.cols());
      for (int i = 0; i < matD2R1R2.rows(); i++) {
        for (int j = 0; j < matD2R1R2.cols(); j++) {
          D2R1R2.set(i, j, matD2R1R2(i, j));
        }
      }
      // Test if the constitutive relations of the resistive elements can be
      // written in input - output form Expensive test - ignored
      // Except the part where matD2R1R2 is updated
      if (posR1.size() > 0 && posR2.size() > 0) {
        auto zeroM = SymbolicMatrix(posR1.size(), posR2.size());
        zeroM.fill(symZero);
        matD2R1R2.block(0, matD2R1R2.cols() - posR2.size() - 1, posR1.size(),
                        posR2.size()) = zeroM;
      }

      SymEngine::DenseMatrix matD2R1R2inv(matD2R1R2.rows(), matD2R1R2.cols());
      D2R1R2.inv(matD2R1R2inv);

      SymbolicMatrix matR2;
      if (posR1.size() == 0) {
        matR2 = matD2R1R2;
      } else {
        // if posR2.size == 0, this is equal to inverse
        matR2 = SymbolicMatrix(matD2R1R2.rows(), matD2R1R2.cols());
        for (int i = 0; i < matD2R1R2.rows(); i++) {
          for (int j = 0; j < matD2R1R2.cols(); j++) {
            matR2(i, j) = matD2R1R2inv.get(i, j);
          }
        }
        if (posR2.size() != 0) {
          SymEngine::DenseMatrix p1(posR1.size(), posR1.size());
          SymEngine::DenseMatrix p1inv(posR1.size(), posR1.size());
          for (int i = 0; i < posR1.size(); i++) {
            for (int j = 0; j < posR1.size(); j++) {
              p1.set(i, j, matD2R1R2(i, j));
            }
          }
          p1.inv(p1inv);
          SymbolicMatrix p1inv_(posR1.size(), posR1.size());
          for (int i = 0; i < posR1.size(); i++) {
            for (int j = 0; j < posR1.size(); j++) {
              p1inv_(i, j) = p1inv.get(i, j);
            }
          }
          SymbolicMatrix p2(posR2.size(), posR2.size());
          int ix = 0;
          int ij = 0;
          for (int i = matD2R1R2.rows() - posR2.size() - 1;
               i < matD2R1R2.rows(); i++) {
            ij = 0;
            for (int j = matD2R1R2.cols() - posR2.size() - 1;
                 j < matD2R1R2.rows(); j++) {
              p2(ix, ij) = matD2R1R2(i, j);
              ij++;
            }
            ix++;
          }
          // ArrayFlatten[{{Inverse[matD2R1R2[[1 ;; Length[posR1], 1 ;;
          // Length[posR1]]]], 0}, {0, matD2R1R2[[-Length[posR2] ;; -1,
          // -Length[posR2] ;; -1]]}}]

          matR2 = SymbolicMatrix(posR1.size() + posR2.size(),
                                 posR1.size() + posR2.size());
          matR2.fill(0);
          matR2.block(0, 0, posR1.size(), posR1.size()) = p1inv_;
          matR2.block(posR1.size(), posR1.size(), posR2.size(), posR2.size()) =
              p2;
        }
      }
      SymbolicMatrix inumR(numR, numR);
      for (int i = 0; i < numR; i++)
        inumR(i, i) = symOne;

      auto matR2zRR = inumR + matR2 * zRR;
      SymEngine::DenseMatrix symMatR2zRR(matR2zRR.rows(), matR2zRR.cols());
      SymEngine::DenseMatrix symMatK(matR2zRR.rows(), matR2zRR.cols());
      for (int i = 0; i < matR2zRR.rows(); i++) {
        for (int j = 0; j < matR2zRR.cols(); j++) {
          symMatR2zRR.set(i, j, matR2zRR(i, j));
        }
      }
      symMatR2zRR.inv(symMatK);
      matK = SymbolicMatrix(matR2zRR.rows(), matR2zRR.cols());
      for (int i = 0; i < matR2zRR.rows(); i++) {
        for (int j = 0; j < matR2zRR.cols(); j++) {
          matK(i, j) = SymEngine::simplify(symMatK.get(i, j));
        }
      }
      matA = matK * matR2 - matR2 * matK.transpose();
      matB = matK * matR2 + matR2 * matK.transpose();
    }
    // case 1 : there are storages, resistors and sources present*
    if (numC > 0 && numR > 0 and numS > 0) {
      matJ = -zCC - SymEngine::parse("1/2") * zCR * matA * zCR.transpose();
      matR = SymEngine::parse("1/2") * zCR * matB * zCR.transpose();
      matG = zCP + SymEngine::parse("1/2") * zCR * matA * zRP;
      matP = SymEngine::parse("-1/2") * zCR * matB * zRP;
      matPT = matP.transpose();
      matM = zPP + SymEngine::parse("1/2") * zRP.transpose() * matA * zRP;
      matS = SymEngine::parse("1/2") * zRP.transpose() * matB * zRP;

      // DHDX = D(funH,{vecX})
    }
    // case 2 : there are resistors and sources but no storages present
    if (numC == 0 && numR > 0 and numS > 0) {
      matM = zPP + SymEngine::parse("1/2") * zRP.transpose() * matA * zRP;
      matS = SymEngine::parse("1/2") * zRP.transpose() * matB * zRP;
    }
    // case 3 : there are storages and sources but no resistors present
    if (numC > 0 && numR == 0 and numS > 0) {
      matJ = -zCC;
      matG = zCP;
      matM = zPP;
      // DHDX = D(funH,{vecX})
    }
    // case 4 : there are sources but no storages and no resistors present
    if (numC == 0 && numR == 0 and numS > 0) {
      matM = zPP;
    }
    // case 5 : there are storages and resistors but no sources present
    if (numC > 0 && numR > 0 and numS == 0) {
      matJ = -zCC - SymEngine::parse("1/2") * zCR * matA * zCR.transpose();
      matR = SymEngine::parse("1/2") * zCR * matB * zCR.transpose();
      // DHDX = D(funH,{vecX})
    }
    // case 6 : there are storages but no resistors and no sources present
    if (numC > 0 && numR == 0 && numS == 0) {
      matJ = -zCC;
      // DHDX = D(funH,{vecX})
    }
    if (numS > 0) {
      ind.clear();
      for (int i = 0; i < numS; i++) {
        ind.push_back(diracEx[1].rows() - numS - 1 + i);
      }
      vecY = diracEx[1](ind, Eigen::placeholders::all);
      vecU = diracEx[3](ind, Eigen::placeholders::all);
    }

    // For computing the Hamiltonian, capacitors have q^2/C, concentrations have
    // RT q (ln (kq) - 1)
    // This is computed using the observation dW = integ[ e dq] limits 0 - q
    // The workdone or change in potential energy due to change in state
    // For capacitors W = Vq; dW = Vdq; dW = (q/C)dq; integral between 0-q ->
    // q^2/2C
#ifdef DEBUG_PHS
    std::cout << "Explicit dirac structure f = D e" << std::endl;
    std::cout << diracEx[1] << std::endl << diracEx[2] << std::endl;
    std::cout << diracEx[3] << std::endl;

    std::cout << "zCC \n" << zCC << std::endl << std::endl;
    std::cout << "zCR \n" << zCR << std::endl << std::endl;
    std::cout << "zCP \n" << zCP << std::endl << std::endl;
    std::cout << "zRR \n" << zRR << std::endl << std::endl;
    std::cout << "zRP \n" << zRP << std::endl << std::endl;
    std::cout << "zPP \n" << zPP << std::endl << std::endl;
    std::cout << "matA \n" << matA << std::endl << std::endl;
    std::cout << "matB \n" << matB << std::endl << std::endl;

    std::cout << "matJ \n" << matJ << std::endl << std::endl;
    std::cout << "matR \n" << matR << std::endl << std::endl;
    std::cout << "matG \n" << matG << std::endl << std::endl;
    std::cout << "matP \n" << matP << std::endl << std::endl;
    std::cout << "matPT \n" << matPT << std::endl << std::endl;
    std::cout << "matM \n" << matM << std::endl << std::endl;
    std::cout << "matS \n" << matS << std::endl << std::endl;
#endif

    result["matA"] = to_json(matA);
    result["matK"] = to_json(matK);
    result["matB"] = to_json(matB);
    result["matJ"] = to_json(matJ);
    result["matR"] = to_json(matR);
    result["matG"] = to_json(matG);
    result["matP"] = to_json(matP);
    result["matM"] = to_json(matM);
    result["matS"] = to_json(matS);
    result["vecY"] = to_json(vecY);
    result["vecU"] = to_json(vecU);
    result["zCC"] = to_json(zCC);
    result["zCR"] = to_json(zCR);
    result["zCP"] = to_json(zCP);
    result["zRR"] = to_json(zRR);
    result["zRP"] = to_json(zRP);
    result["zPP"] = to_json(zPP);

    nlohmann::json ds;
    ds["lhs"] = to_json(diracEx[1]);
    ds["ds"] = to_json(diracEx[2]);
    ds["rhs"] = to_json(diracEx[3]);
    result["ExplicitDS"] = ds;

    // Port hamiltonian
    SymbolicMatrix matJR = matJ - matR;
    SymbolicMatrix matGP = matG - matP;
    SymbolicMatrix matGt = matG.transpose();
    SymbolicMatrix matGtPt = matGt - matPT;
    SymbolicMatrix matMS = matM + matS;
    nlohmann::json phs;
    phs["matJR"] = to_json(matJR);
    phs["matGP"] = to_json(matGP);
    phs["matGt"] = to_json(matGt);
    phs["matGtPt"] = to_json(matGtPt);
    phs["matMS"] = to_json(matMS);
    result["portHamiltonianMatrices"] = phs;
    result["success"] = true;
  } catch (BGException &ex) {
    result["error"] = ex.what();
  }
  if (warnings.size())
    result["warnings"] = warnings;
  std::ofstream file("d:/Temp/LatexRenderer/phs.json");
  file << result;
  return result;
}