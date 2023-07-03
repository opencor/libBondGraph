

nlohmann::json to_json(const RCPLIB::RCP<const SymEngine::DenseMatrix> &mat) {
  nlohmann::json js;
  js["rows"] = mat->nrows();
  js["cols"] = mat->ncols();
  std::vector<std::string> elems;
  for (int i = 0; i < mat->nrows(); i++) {
    for (int j = 0; j < mat->ncols(); j++) {
      SymEngine::Expression simpE = SymEngine::simplify(mat->get(i, j));
      // elems.push_back(SymEngine::mathml(simpE));
      elems.push_back(SymEngine::latex(simpE));
    }
  }
  js["elements"] = elems;
  return js;
}

nlohmann::json to_json(SymEngine::DenseMatrix &mat) {
  nlohmann::json js;
  js["rows"] = mat.nrows();
  js["cols"] = mat.ncols();
  std::vector<std::string> elems;
  for (int i = 0; i < mat.nrows(); i++) {
    for (int j = 0; j < mat.ncols(); j++) {
      SymEngine::Expression simpE = SymEngine::simplify(mat.get(i, j));
      // elems.push_back(SymEngine::mathml(simpE));
      elems.push_back(SymEngine::latex(simpE));
    }
  }
  js["elements"] = elems;
  return js;
}

nlohmann::json to_json(const SymEngine::vec_basic &mat) {
  nlohmann::json js;
  js["rows"] = mat.size();
  js["cols"] = 1;
  std::vector<std::string> elems;
  for (int i = 0; i < mat.size(); i++) {
    SymEngine::Expression simpE = SymEngine::simplify(mat[i]);
    // elems.push_back(SymEngine::mathml(simpE));
    elems.push_back(SymEngine::latex(simpE));
  }
  js["elements"] = elems;
  return js;
}

nlohmann::json BondGraph::computePortHamiltonian() {
  nlohmann::json result; // Compose results as we progress
  result["success"] = false;
  try {
    SymEngine::vec_basic dofs;
    // dx,e,f,x,u
    SymEngine::vec_basic dx;
    SymEngine::vec_basic e;
    SymEngine::vec_basic f;
    SymEngine::vec_basic x;
    SymEngine::vec_basic u;
    std::unordered_map<std::string, std::string> nameMap;
    std::unordered_map<std::string, RCPLIB::RCP<const SymEngine::Basic>>
        nonlinearTermSubs;
    std::unordered_map<std::string, std::tuple<std::string, std::string, char>>
        dimensions; // variable name, value, state, control or parameter
    std::unordered_map<std::string, std::vector<nlohmann::json>> annotations;
    std::unordered_map<size_t, SymEngine::map_basic_basic>
        portHamiltonianCoordinates;

    std::vector<std::string> hamiltonianString;
    std::vector<std::string> stateVariables;
    SymEngine::vec_basic parameterDivisor; // Parameter by which JR matrix entry
                                           // should be divided by as Derivative
                                           // of Hamiltonian will have it
    std::vector<bool> chemicalStorage;
    // Global dof matrix
    // Count the number of element in each group
    mUcount = 0; // Number of source bonds
    mScount = 0; // Number of storage bonds with integral causality
    mRcount = 0; // Number of dissipative bonds
    mJcount = 0; // Number of junction bonds
    mSources.clear();
    //
    // Find all the connected components

    std::vector<RCPLIB::RCP<BGElement>> connectedComponents;
    std::map<std::string, RCPLIB::RCP<BGElement>> compIDMap;
    for (auto &bd : mBonds) {
      auto from = bd->getFromPort()->getComponent();
      auto to = bd->getToPort()->getComponent();
      compIDMap[from->getId()] = from;
      compIDMap[to->getId()] = to;
    }
    for (auto &c : mComponents) { // Maintain order
      auto id = c->getId();
      if (compIDMap.find(id) != compIDMap.end()) {
        connectedComponents.push_back(compIDMap[id]);
      }
    }

    logDebug(" Connected components ", connectedComponents.size(),
             " num bonds ", mBonds.size());
    long int dofID = 0;
    int portID = 0;
    std::ostringstream ss, dss, ess, fss, pss;

    std::unordered_map<std::string, unsigned int> uniqueElementNames;

    // Do storage, dissipative, junction, source
    // Assign unique dof values to the components to be pulled into a system of
    // equations
    for (auto &mc_ : connectedComponents) {
      const RCPLIB::RCP<BondGraphElementBase> mc =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(mc_);

      if (mc->getComponentGroup() == eS) {
        // std::string eName = mc->getName();
        std::string eName = mc->getVariableName();
        if (uniqueElementNames.find(eName) != uniqueElementNames.end()) {
          unsigned int ec = uniqueElementNames[eName];
          uniqueElementNames[eName] = ec + 1;
          eName = eName + "_" + std::to_string(ec);
        } else {
          uniqueElementNames[eName] = 1;
        }

        auto values = mc->values();
        ss.str("");
        ss.clear();
        dss.str("");
        dss.clear();
        ess.str("");
        ess.clear();
        fss.str("");
        fss.clear();
        std::string stateName = std::get<0>(values[0]);
        auto ploc = stateName.rfind("_");
        if (ploc != std::string::npos) {
          ss << stateName.substr(0, ploc + 1);
        } else {
          ss << stateName;
        }
        ss << dofID;
        mc->setDof(dofID);
        x.push_back(SymEngine::symbol(ss.str()));
        dss << "dot_" << ss.str();
        dx.push_back(SymEngine::symbol(dss.str()));
        nameMap[ss.str()] =
            stateName.substr(0, ploc + 1) + "of_" + eName; //+"_"+ss.str();
        nameMap[dss.str()] = "dot_" + stateName.substr(0, ploc + 1) + "of_" +
                             eName; //+"_"+ss.str();//dss.str();

        ess << "e_" << portID;
        fss << "f_" << portID;
        e.push_back(SymEngine::symbol(ess.str()));
        f.push_back(SymEngine::symbol(fss.str()));
        mc->getPorts(0)->setDofIndex(portID);

        auto un = std::get<1>(values[0])->units;
        auto vn = std::get<1>(values[0])->value;
        // auto mult = un.multiplier();
        // auto baseU = un.base_units();
        nameMap[ess.str()] = eName + "_" + ess.str();
        nameMap[fss.str()] = eName + "_" + fss.str();

        // Create hamiltonian
        // For computing the Hamiltonian, capacitors have q^2/C, concentrations
        // have RT q (ln (kq) - 1) This is computed using the observation dW =
        // integ[ e dq] limits 0 - q The workdone or change in potential energy
        // due to change in state For capacitors W = Vq; dW = Vdq; dW = (q/C)dq;
        // integral between 0-q -> q^2/2C

        switch (mc->getType()) {
        case eCapacitance:
        case eInductance: {
          int numStates = mc->getNumStates();
          auto values = mc->values();
          std::string state =
              std::get<1>(values[numStates - 1])->name + std::to_string(dofID);
          stateVariables.push_back(state);
          std::string parm = std::get<1>(values[numStates])->name;
          std::string energy = "(1/2)*(1/" + parm + ")*" + state + "**2";
          hamiltonianString.push_back(energy);
          parameterDivisor.push_back(SymEngine::parse(parm));
          chemicalStorage.push_back(false);
          break;
        }
        case bConcentration: {
          int numStates = mc->getNumStates();
          auto values = mc->values();
          std::string state =
              std::get<1>(values[numStates - 1])->name + std::to_string(dofID);
          stateVariables.push_back(state);
          std::string parm = std::get<1>(values[numStates + 2])->name;
          std::string pR, pT;
          auto Rv = std::get<1>(values[numStates]);
          if (Rv->universalConstant) {
            pR = Rv->name;
          } else {
            // This needs to be checked as param names are updated
            pR = Rv->name + "_" + std::to_string(dofID);
          }
          auto Rt = std::get<1>(values[numStates + 1]);
          if (Rt->universalConstant) {
            pT = Rt->name;
          } else {
            // This needs to be checked as param names are updated
            pT = Rt->name + "_" + std::to_string(dofID);
          }
          // RTq_0(ln(kq_0)−1)
          std::string energy = pR + "*" + pT + "*" + state + "*(log(" + parm +
                               "*" + state + ")-1.0)";
          parameterDivisor.push_back(SymEngine::parse(parm));
          chemicalStorage.push_back(true);
          hamiltonianString.push_back(energy);
          break;
        }
        default: {
          parameterDivisor.push_back(SymEngine::one);
        }
        }

        ++portID;
        ++dofID;
        ++mScount;

        for (int i = mc->getNumStates(); i < values.size(); i++) {
          // Use the actual name instead of prefix
          std::string pname = std::get<1>(values[i])->name;
          if (!std::get<1>(values[i])->universalConstant) {
            // std::cout<<pname<<" ->
            // "<<std::get<1>(values[i])->prefix+"_of_"+eName<<std::endl;
            nameMap[pname] = std::get<1>(values[i])->prefix + "_of_" + eName;
          }
        }
      }
      if (mc->getComponentGroup() == ePH) {
        // std::string eName = mc->getName();
        logWarn("Current implementation does not implement PHS");
        std::string eName = mc->getVariableName();
        if (uniqueElementNames.find(eName) != uniqueElementNames.end()) {
          unsigned int ec = uniqueElementNames[eName];
          uniqueElementNames[eName] = ec + 1;
          eName = eName + "_" + std::to_string(ec);
        } else {
          uniqueElementNames[eName] = 1;
        }
        auto values = mc->values();
        SymEngine::map_basic_basic globalCoordinates;
        // The entire ph has one dofID
        ++dofID;
        mc->setDof(dofID);
        for (int st = 0; st < mc->getNumStates(); st++) {
          ss.str("");
          ss.clear();
          dss.str("");
          dss.clear();
          ess.str("");
          ess.clear();
          fss.str("");
          fss.clear();
          std::string stateName = std::get<0>(values[st]);
          auto ploc = stateName.rfind("_");
          if (ploc != std::string::npos) {
            ss << stateName.substr(0, ploc + 1);
          } else {
            ss << stateName;
          }
          ss << dofID;

          x.push_back(SymEngine::symbol(ss.str()));
          dss << "dot_" << ss.str();
          dx.push_back(SymEngine::symbol(dss.str()));
          nameMap[ss.str()] =
              stateName.substr(0, ploc + 1) + "of_" + eName; //+"_"+ss.str();
          nameMap[dss.str()] = "dot_" + stateName.substr(0, ploc + 1) + "of_" +
                               eName; //+"_"+ss.str();

          ess << "e_" << portID;
          fss << "f_" << portID;
          e.push_back(SymEngine::symbol(ess.str()));
          f.push_back(SymEngine::symbol(fss.str()));
          mc->getPorts(st)->setDofIndex(portID);
          nameMap[ess.str()] = eName + "_" + ess.str();
          nameMap[fss.str()] = eName + "_" + fss.str();

          auto un = std::get<1>(values[st])->units;
          auto vn = std::get<1>(values[st])->value;

          globalCoordinates[SymEngine::parse(stateName)] =
              SymEngine::parse(ss.str());
          globalCoordinates[SymEngine::parse("dot_" + stateName)] =
              SymEngine::parse("dot_" + ss.str());
          globalCoordinates[SymEngine::parse("e_" + std::to_string(st))] =
              SymEngine::parse(ess.str());
          globalCoordinates[SymEngine::parse("f_" + std::to_string(st))] =
              SymEngine::parse(fss.str());

          ++portID;
          ++mScount;
        }
        // Dimensions for the parameters
        for (int i = mc->getNumStates(); i < values.size(); i++) {
          // Use the actual name instead of prefix
          std::string pname = std::get<1>(values[i])->name;
          globalCoordinates[SymEngine::parse(std::get<0>(values[i]))] =
              SymEngine::parse(pname);
          if (!std::get<1>(values[i])->universalConstant) {
            nameMap[pname] = std::get<1>(values[i])->prefix + "_of_" + eName;
          }
        }
        portHamiltonianCoordinates[dofID] = globalCoordinates;
      }
    }
    // Do for resistors
    for (auto &mc_ : connectedComponents) {
      const RCPLIB::RCP<BondGraphElementBase> mc =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(mc_);
      if (mc->getComponentGroup() == eR) {
        // std::string eName = mc->getName();
        std::string eName = mc->getVariableName();
        if (uniqueElementNames.find(eName) != uniqueElementNames.end()) {
          unsigned int ec = uniqueElementNames[eName];
          uniqueElementNames[eName] = ec + 1;
          eName = eName + "_" + std::to_string(ec);
        } else {
          uniqueElementNames[eName] = 1;
        }
        auto ports = mc->getPorts();
        // Handle reactions
        for (int pi = 0; pi < ports.size(); pi++) {
          ess.str("");
          ess.clear();
          fss.str("");
          fss.clear();
          ess << "e_" << portID;
          fss << "f_" << portID;
          e.push_back(SymEngine::symbol(ess.str()));
          f.push_back(SymEngine::symbol(fss.str()));
          nameMap[ess.str()] = eName + "_" + ess.str();
          nameMap[fss.str()] = eName + "_" + fss.str();
          ports[pi]->setDofIndex(portID);
          ++portID;
        }
        mc->setDof(dofID);
        ++dofID;
        ++mRcount;

        auto values = mc->values();
        for (int i = mc->getNumStates(); i < values.size(); i++) {
          // Use the actual name instead of prefix
          std::string pname = std::get<1>(values[i])->name;
          if (!std::get<1>(values[i])->universalConstant) {
            nameMap[pname] = std::get<1>(values[i])->prefix + "_of_" + eName;
          }
        }
      }
    }
    // Do for sources
    for (auto &mc_ : connectedComponents) {
      const RCPLIB::RCP<BondGraphElementBase> mc =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(mc_);
      // std::string eName = mc->getName();
      if (mc->getComponentGroup() == eU) {
        std::string eName = mc->getVariableName();
        if (uniqueElementNames.find(eName) != uniqueElementNames.end()) {
          unsigned int ec = uniqueElementNames[eName];
          uniqueElementNames[eName] = ec + 1;
          eName = eName + "_" + std::to_string(ec);
        } else {
          uniqueElementNames[eName] = 1;
        }
        ess.str("");
        ess.clear();
        fss.str("");
        fss.clear();
        ess << "e_" << portID;
        fss << "f_" << portID;
        e.push_back(SymEngine::symbol(ess.str()));
        f.push_back(SymEngine::symbol(fss.str()));
        mc->setDof(dofID);
        mc->getPorts(0)->setDofIndex(portID);

        nameMap[ess.str()] = eName + "_" + ess.str();
        nameMap[fss.str()] = eName + "_" + fss.str();

        auto mcValues = mc->values();
        for (int pi = mc->getNumStates(); pi < mcValues.size(); pi++) {
          ess.str("");
          ess.clear();
          ess << std::get<0>(mcValues[pi]) << "_" << dofID;
          // ess << std::get<0>(mcValues[pi]) << "_of_" << eName;
          u.push_back(SymEngine::symbol(ess.str()));

          auto un = std::get<1>(mcValues[pi])->units;
          auto vn = std::get<1>(mcValues[pi])->value;

          nameMap[ess.str()] = std::get<0>(mcValues[pi]) + "_of_" + eName;
        }
        mSources.push_back(mc);
        ++portID;
        ++dofID;
        ++mUcount;
      }
    }
    // For junctions
    for (auto &mc_ : connectedComponents) {
      const RCPLIB::RCP<BondGraphElementBase> mc =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(mc_);
      if (mc->getComponentGroup() == eJ) {
        // std::string eName = mc->getName();
        std::string eName = mc->getVariableName();
        if (uniqueElementNames.find(eName) != uniqueElementNames.end()) {
          unsigned int ec = uniqueElementNames[eName];
          uniqueElementNames[eName] = ec + 1;
          eName = eName + "_" + std::to_string(ec);
        } else {
          uniqueElementNames[eName] = 1;
        }
        auto ports = mc->getPorts();
        for (int i = 0; i < ports.size(); i++) {
          ess.str("");
          ess.clear();
          fss.str("");
          fss.clear();
          ess << "e_" << portID;
          fss << "f_" << portID;
          e.push_back(SymEngine::symbol(ess.str()));
          f.push_back(SymEngine::symbol(fss.str()));
          nameMap[ess.str()] = eName + "_" + ess.str();
          nameMap[fss.str()] = eName + "_" + fss.str();

          ports[i]->setDofIndex(portID);
          ++portID;
        }
        mc->setDof(dofID);
        ++dofID;
        ++mJcount;
      }
    }

    dofs.insert(dofs.end(), dx.begin(), dx.end());
    for (int ec = 0; ec < e.size(); ec++) {
      dofs.push_back(e[ec]);
      dofs.push_back(f[ec]);
    }
    dofs.insert(dofs.end(), x.begin(), x.end());
    dofs.insert(dofs.end(), u.begin(), u.end());
    unsigned int numStates = x.size();
    unsigned int numBonds = e.size();
    unsigned int numControlVariables = u.size();
    auto coordinates = SymEngine::DenseMatrix(dofs);

    // JS matrix
    unsigned int bc = 0;
    auto offset = mScount;
    std::vector<int> ix, iy, dir;
    for (auto &bd : mBonds) {
      int j = bd->getFromPort()->dofIndex();
      int k = bd->getToPort()->dofIndex();
      ix.push_back(2 * bc);
      iy.push_back(j * 2 + offset);
      dir.push_back(-1);
      ix.push_back(2 * bc);
      iy.push_back(k * 2 + offset);
      dir.push_back(1);
      ix.push_back(2 * bc + 1);
      iy.push_back(j * 2 + offset + 1);
      dir.push_back(1);
      ix.push_back(2 * bc + 1);
      iy.push_back(k * 2 + offset + 1);
      dir.push_back(1);
      bc++;
    }

    unsigned int rows = 2 * bc;
    // Handle constitutive relations
    // Create coordinate map
    coordinateMap.clear();

    long int di = 0;
    for (auto &c : dofs) {
      ss.str("");
      ss.clear();
      ss << *c;
      coordinateMap[ss.str()] = di++;
    }
    std::unordered_map<unsigned int, std::vector<ExpressionTermsMap>> cIndexes;
    for (auto &mc_ : connectedComponents) {
      const RCPLIB::RCP<BondGraphElementBase> mc =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(mc_);
      std::vector<ExpressionTermsMap> cx =
          getLinearAndNonlinearTerms(mc, coordinateMap);
      cIndexes[rows] = cx;

      for (int rc = 0; rc < mc->getConstitutiveEquations().size(); rc++) {
        mc->setConstitutiveEqIndex(rows, rc);
        rows++;
      }
    }

    SymEngine::DenseMatrix linOp(rows, dofs.size());
    SymEngine::DenseMatrix nonlinearTerms(rows, 1);
    zeros(nonlinearTerms); // Initialize the matices with zeros
    zeros(linOp);
    for (int i = 0; i < ix.size(); i++) {
      linOp.set(ix[i], iy[i], SymEngine::integer(dir[i]));
    }

    // Constitutive equations
    for (auto &mp : cIndexes) {
      auto rowIx = mp.first;
      for (auto &mx : mp.second) {
        for (int i = 0; i < mx.indexes.size(); i++) {
          auto col = mx.indexes[i];
          auto cef = mx.coefficients[i];
          if (col > -1) {
            linOp.set(rowIx, col, cef);
          } else {
            nonlinearTerms.set(
                rowIx, 0, SymEngine::add(nonlinearTerms.get(rowIx, 0), cef));
          }
        }
        rowIx++;
      }
    }

    // Optimise findSubstitutions for speed
    auto snf = getSmithNormalForm(linOp, nonlinearTerms);
    // Add constraints from Port Hamiltionians
    for (auto &mc_ : connectedComponents) {
      const RCPLIB::RCP<BondGraphElementBase> mc =
          RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(mc_);
      if (mc->getComponentGroup() == ePH) {
        RCPLIB::RCP<PortHamiltonian> ph =
            RCPLIB::rcp_dynamic_cast<PortHamiltonian>(mc_);
        const std::vector<std::string> &cons = ph->getConstraints();
        for (auto &c : cons) {
          // Go from local to global coordinates
          auto cexp = SymEngine::parse(c);
          long int id = ph->getDof();
          auto &subs = portHamiltonianCoordinates[id];
          snf.constraints.push_back(cexp->subs(subs));
        }
      }
    }

    // Optimise findSubstitutions for speed
    auto subsExprs =
        findSubstitutions(snf.linearOp, snf.nonlinearOp, snf.constraints,
                          coordinates, numStates, numBonds);
    substituteValues(snf.nonlinearOp, subsExprs);
    substituteValues(snf.constraints, subsExprs);

    // Process constraints
    auto sn = process_constraints(snf, coordinates, numStates, numBonds,
                                  numControlVariables);
    auto subsExprsC =
        findSubstitutions(sn.linearOp, sn.nonlinearOp, sn.constraints,
                          coordinates, numStates, numBonds);
    substituteValues(sn.nonlinearOp, subsExprsC);
    substituteValues(sn.constraints, subsExprsC);

    SymEngine::DenseMatrix &linearOp = sn.linearOp;
    SymEngine::DenseMatrix &nonlinearOp = sn.nonlinearOp;
    SymEngine::vec_basic &constraints = sn.constraints;

    int rows_added = 0;
    std::vector<int> added_cvs;
    SymEngine::map_basic_basic cv_diff_dict;
    std::vector<std::tuple<int, int, SymEngine::RCP<const SymEngine::Basic>>>
        lin_dict;

    // Create coordinate map
    coordinateMap.clear();
    di = 0;
    for (auto &c : dofs) {
      ss.str("");
      ss.clear();
      ss << *c;
      std::string dofname = ss.str(); // Eliminate memory loss issue
      coordinateMap[dofname] = di++;
    }

    bool solvable = true; // Flag to check if the equations are solvable.. Some
                          // incorrect bg formulations will lead to state and
                          // bond variables being solved to zero
    // - Linear constraints; ie Lx = 0
    // - Nonlinear Constraints Lx + F(x) = 0
    //
    // Linear constraints are rows with more than 1 non-zero
    // that are not in the derivative subspace, and have a zero nonlinear part
    offset = 2 * numBonds + numStates;
    if (linearOp.nrows() < offset) {
      SymEngine::set_basic c_atoms;
      for (int i = 0; i < coordinates.nrows(); i++) {
        auto atoms =
            SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(
                *coordinates.get(i, 0));
        c_atoms.insert(atoms.begin(), atoms.end());
      }
      for (int row = offset - 1; row > linearOp.nrows(); row--) {
        auto atoms =
            SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(
                *nonlinearOp.get(row, 0));
        auto a_c = setIntersection(c_atoms, atoms);
        if (a_c.size() == 0) {
          SymEngine::DenseMatrix linRow(1, linearOp.ncols());
          linearOp.submatrix(linRow, row, 0, row, linearOp.ncols() - 1);
          if (!is_true(linRow.is_zero())) {
            for (int idx = 0; idx < numStates; idx++) {
              auto v = linearOp.get(row, idx + offset);
              if (!SymEngine::eq(*v, *SymEngine::zero)) {
                lin_dict.push_back(std::make_tuple(rows_added, idx, v));
              }
            }
            for (int idx = 0; idx < numControlVariables; idx++) {
              auto v = linearOp.get(row, idx + offset + numStates);
              if (!SymEngine::eq(*v, *SymEngine::zero)) {
                lin_dict.push_back(std::make_tuple(rows_added, idx, v));
              }
            }
          }
        }
      }
    }

    std::vector<std::tuple<int, int, SymEngine::RCP<const SymEngine::Basic>>>
        cv_dict;
    rows = linearOp.nrows(); // Only parse for existing rows, not the ones that
                             // will be added
    for (int row = offset; row < rows; row++) {
      auto nonlinear_constraint = nonlinearOp.get(row, 0);
      // Size change during the call so do not allocate outside
      SymEngine::DenseMatrix linearOpSub(1, linearOp.ncols() - offset);
      linearOp.submatrix(linearOpSub, row, offset, row, linearOp.ncols() - 1);
      if (is_true(linearOpSub.is_zero()) &&
          SymEngine::eq(*nonlinear_constraint, *SymEngine::zero)) {
        continue;
      }
      // Size change during the call so do not allocate outside
      SymEngine::DenseMatrix state_constraint(1, numStates);
      SymEngine::DenseMatrix control_constraint(1, linearOp.ncols() - offset -
                                                       numStates);
      SymEngine::DenseMatrix sj(1, offset + numControlVariables);

      linearOp.submatrix(state_constraint, row, offset, row,
                         offset + numStates - 1);
      linearOp.submatrix(control_constraint, row, offset + numStates, row,
                         linearOp.ncols() - 1);
      zeros(sj);

      state_constraint.row_join(sj);

      cv_dict.clear();
      if (!is_true(control_constraint.is_zero())) {
        for (int cv_col = 0; cv_col < control_constraint.ncols(); cv_col++) {
          auto con = control_constraint.get(0, cv_col);
          if (SymEngine::eq(*con, *SymEngine::zero)) {
            continue;
          }
          int idx = added_cvs.size();
          if (std::find(added_cvs.begin(), added_cvs.end(), cv_col) ==
              added_cvs.end()) {
            added_cvs.push_back(cv_col);
            SymEngine::DenseMatrix loj(linearOp.nrows(), 1);
            zeros(loj);
            linearOp.row_join(loj);
            auto coord = coordinates.get(offset + numStates + cv_col, 0);
            ss.str("");
            ss.clear();
            ss << "dot_" << *coord;
            SymEngine::DenseMatrix coj(1, 1);
            coj.set(0, 0, SymEngine::parse(ss.str()));
            coordinates.col_join(coj);
            numControlVariables++;
          } else {
            idx = std::find(added_cvs.begin(), added_cvs.end(), cv_col) -
                  added_cvs.begin();
          }
          cv_dict.push_back(std::make_tuple(0, idx, con));
        }
      }

      if (added_cvs.size() > 0) {
        SymEngine::DenseMatrix cons(1, added_cvs.size());
        for (auto &t : cv_dict) {
          int i, j;
          SymEngine::RCP<const SymEngine::Basic> v;
          std::tie(i, j, v) = t;
          cons.set(i, j, v);
        }
        state_constraint.row_join(cons);
      }
      SymEngine::vec_basic jac_dx;
      bool secondOrderConstraint = false;
      for (int ci = 0; ci < numStates; ci++) {
        auto c = coordinates.get(ci, 0);
        ss.str("");
        ss.clear();
        ss << *c;
        SymEngine::RCP<const SymEngine::Symbol> csym =
            SymEngine::symbol(ss.str());
        auto ncd = nonlinear_constraint->diff(csym);
        jac_dx.push_back(ncd);
        if (!SymEngine::eq(*ncd, *SymEngine::zero)) {
          secondOrderConstraint = true;
        }
      }
      SymEngine::vec_basic jac_junction;
      bool firstOrderJunctionConstraint = false;
      for (int ci = numStates; ci < offset; ci++) {
        auto c = coordinates.get(ci, 0);
        ss.str("");
        ss.clear();
        ss << *c;
        SymEngine::RCP<const SymEngine::Symbol> csym =
            SymEngine::symbol(ss.str());
        auto ncd = nonlinear_constraint->diff(csym);
        jac_junction.push_back(ncd);
        if (!SymEngine::eq(*ncd, *SymEngine::zero)) {
          firstOrderJunctionConstraint = true;
        }
      }
      SymEngine::vec_basic jac_x;
      bool secondOrderXConstraint = false;
      for (int ci = offset; ci < offset + numStates; ci++) {
        auto c = coordinates.get(ci, 0);
        ss.str("");
        ss.clear();
        ss << *c;
        SymEngine::RCP<const SymEngine::Symbol> csym =
            SymEngine::symbol(ss.str());
        auto ncd = nonlinear_constraint->diff(csym);
        jac_x.push_back(ncd);
        if (!SymEngine::eq(*ncd, *SymEngine::zero)) {
          secondOrderXConstraint = true;
        }
      }
      SymEngine::vec_basic jac_cv;
      bool firstOrderControlConstraint = false;
      for (int ci = offset + numStates; ci < coordinates.nrows(); ci++) {
        auto c = coordinates.get(ci, 0);
        ss.str("");
        ss.clear();
        ss << *c;
        SymEngine::RCP<const SymEngine::Symbol> csym =
            SymEngine::symbol(ss.str());
        auto ncd = nonlinear_constraint->diff(csym);
        jac_cv.push_back(ncd);
        if (!SymEngine::eq(*ncd, *SymEngine::zero)) {
          firstOrderControlConstraint = true;
        }
      }
      auto nlin_row = SymEngine::parse("0");

      if (secondOrderXConstraint) {
        auto sumxy = SymEngine::mul(jac_x[0], coordinates.get(0, 0));
        for (int x = 1; x < numStates; x++) {
          sumxy = SymEngine::add(
              sumxy, SymEngine::mul(jac_x[x], coordinates.get(x, 0)));
        }
        auto expr = SymEngine::expand(sumxy);
        SymEngine::RCP<const SymEngine::Basic> num, denom;
        SymEngine::as_numer_denom(expr, SymEngine::outArg(num),
                                  SymEngine::outArg(denom));
        if (is_true(state_constraint.is_zero())) {
          std::map<long int, SymEngine::RCP<const SymEngine::Basic>> ld;
          SymEngine::RCP<const SymEngine::Basic> nl;
          SymEngine::map_basic_basic dummy;
          ss.str("");
          ss.clear();
          ss << *num;
          std::string numer = ss.str();
          std::tie(ld, nl) = getLinearCoefficientsAndNonlinearTerms(
              numer, coordinateMap, dummy);
          for (auto &kv : ld) {
            auto ov = state_constraint.get(0, kv.first);
            state_constraint.set(0, kv.first, SymEngine::add(ov, kv.second));
          }
          nlin_row = SymEngine::add(nlin_row, nl);
        } else {
          nlin_row = SymEngine::add(nlin_row, expr);
        }
      }
      SymEngine::DenseMatrix nadd(1, 1);
      nadd.set(0, 0, nlin_row);
      nonlinearOp.col_join(nadd);
      linearOp.col_join(state_constraint);
      rows_added++;
    }

    if (rows_added > 0) {
      auto snf = getSmithNormalForm(linearOp, nonlinearOp);
      linearOp = snf.linearOp;
      nonlinearOp = snf.nonlinearOp;
      constraints = snf.constraints;
    }

    // Simplify Exp(Log(x)) terms in nonlinear Op
    for (int nr = 0; nr < nonlinearOp.nrows(); nr++) {
      if (!SymEngine::eq(*nonlinearOp.get(nr, 0), *SymEngine::zero)) {
        nonlinearOp.set(nr, 0,
                        SymEngine::simplifyExpLog(nonlinearOp.get(nr, 0)));
      }
    }

    auto System = SymEngine::DenseMatrix(linearOp.nrows(), coordinates.ncols());
    mul_dense_dense(linearOp, coordinates, System);

    std::ostringstream hms;
    for (auto &e : hamiltonianString) {
      hms << e + " + ";
    }
    std::string hmString = hms.str();
    auto hmexp = SymEngine::simplify(SymEngine::expand(
        SymEngine::parse(hmString.substr(0, hmString.size() - 3))));
    // Derivatives
    SymEngine::DenseMatrix derivatives(stateVariables.size(), 1);
    const RCPLIB::RCP<const SymEngine::Basic> &ham = hmexp;
    for (size_t svi = 0; svi < stateVariables.size(); svi++) {
      RCPLIB::RCP<const SymEngine::Symbol> sv =
          SymEngine::symbol(stateVariables[svi]);
      derivatives.set(svi, 0, SymEngine::diff(ham, sv, true));
    }

    // For linear terms divide them by the parameter
    offset = 2 * numBonds + numStates;
    for (int xs = 0; xs < numStates; xs++) {
      for (size_t xs1 = 0; xs1 < numStates; xs1++) {
        auto lhs = linearOp.get(xs, offset + xs1);
        if (!SymEngine::eq(*lhs, *SymEngine::zero)) {
          SymEngine::map_basic_basic subt;
          subt[parameterDivisor[xs1]] = SymEngine::one;
          linearOp.set(xs, offset + xs1, SymEngine::expand(lhs->subs(subt)));
        }
      }
    }
    // Loop through nonlinear terms for the derivatives
    // Assign terms to the linearop matrix at the state col

    for (int i = 0; i < numStates; i++) {
      auto rhs = nonlinearOp.get(i, 0);
      if (!SymEngine::eq(*rhs, *SymEngine::zero)) {
        if (SymEngine::is_a<SymEngine::Add>(*rhs)) {
          auto terms = rhs->get_args();
          for (const auto &t : terms) {
            auto atoms =
                SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(
                    *t);
            if (atoms.find(x[i]) != atoms.end()) {
              SymEngine::map_basic_basic subt;
              subt[x[i]] = SymEngine::one;
              subt[parameterDivisor[i]] = SymEngine::one;
              auto nexp = SymEngine::expand(t->subs(subt));
              auto existing = linearOp.get(i, offset + i);
              SymEngine::vec_basic aterms;
              aterms.push_back(existing);
              ss.str("");
              ss.clear();
              ss << *nexp << "*ExpRT";
              aterms.push_back(SymEngine::parse(ss.str()));
              linearOp.set(i, offset + i, SymEngine::add(aterms));
            }
          }
        }
      }
    }

    SymEngine::DenseMatrix JR(numStates, numStates);
    SymEngine::DenseMatrix E(numStates, numStates);
    SymEngine::eye(E);
    for (size_t xs = 0; xs < numStates; xs++) {
      for (size_t xs1 = 0; xs1 < numStates; xs1++) {
        JR.set(xs, xs1, linearOp.get(xs, xs1 + offset));
      }
    }

    // Extract J, R and Q matrix
    SymEngine::DenseMatrix JR2(numStates, numStates);
    SymEngine::DenseMatrix J(numStates, numStates);
    SymEngine::DenseMatrix R(numStates, numStates);
    SymEngine::DenseMatrix Q(numStates, numStates);
    SymEngine::eye(Q);
    for (int ci = 0; ci < numStates; ci++) {
      if (!chemicalStorage[ci]) {
        Q.set(ci, ci, SymEngine::div(SymEngine::one, parameterDivisor[ci]));
      }
    }

    SymEngine::DenseMatrix JRT(numStates, numStates);
    JR.transpose(JRT);
    // Symmetric part is (A + A^T)/2
    JR.add_matrix(JRT, JR2);
    JR2.mul_scalar(SymEngine::parse("1/2"), R);
    // Skew symmetric  part is (A - A^T)/2
    JRT.mul_scalar(SymEngine::parse("-1"), JR2);
    JR.add_matrix(JR2, JRT);
    JRT.mul_scalar(SymEngine::parse("1/2"), J);

    size_t ucols = mSources.size() == 0 ? 1 : mSources.size();
    SymEngine::DenseMatrix B(numStates, ucols);
    SymEngine::DenseMatrix U(1, numStates);
    SymEngine::zeros(B);
    SymEngine::zeros(U);
    for (size_t xs = 0; xs < numStates; xs++) {
      for (size_t bs = 0; bs < mSources.size(); bs++) {
        auto entry = linearOp.get(xs, bs + offset + numStates);
        B.set(xs, bs, entry);
        if (!SymEngine::eq(*entry, *SymEngine::zero)) {
          U.set(0, xs, coordinates.get(bs + offset + numStates, 0));
        }
      }
    }

    // Do name mapping
    SymEngine::map_basic_basic nameMapSubs;
    for (auto &c : nameMap) {
      nameMapSubs[SymEngine::symbol(c.first)] = SymEngine::symbol(c.second);
    }

    for (int i = 0; i < J.nrows(); i++) {
      for (int j = 0; j < J.ncols(); j++) {
        J.set(i, j, J.get(i, j)->subs(nameMapSubs));
      }
    }

    for (int i = 0; i < R.nrows(); i++) {
      for (int j = 0; j < R.ncols(); j++) {
        R.set(i, j, R.get(i, j)->subs(nameMapSubs));
      }
    }

    for (int i = 0; i < U.nrows(); i++) {
      for (int j = 0; j < U.ncols(); j++) {
        U.set(i, j, U.get(i, j)->subs(nameMapSubs));
      }
    }

    for (int i = 0; i < derivatives.nrows(); i++) {
      for (int j = 0; j < derivatives.ncols(); j++) {
        derivatives.set(i, j, derivatives.get(i, j)->subs(nameMapSubs));
      }
    }

    for (int i = 0; i < x.size(); i++) {
      x[i] = x[i]->subs(nameMapSubs);
    }

    hmexp->subs(nameMapSubs);

    // Completed NameMap

    nlohmann::json phs;
    phs["matJ"] = to_json(J);
    phs["matR"] = to_json(R);
    phs["matB"] = to_json(B);
    phs["u"] = to_json(U);
    phs["matE"] = to_json(E);
    phs["matQ"] = to_json(Q);

    std::cout << *hmexp << std::endl;
    std::cout << "J" << std::endl;
    std::cout << J << std::endl;
    std::cout << "R" << std::endl;
    std::cout << R << std::endl;

    result["stateVector"] = to_json(x);
    result["Hderivatives"] = to_json(derivatives);
    result["hamiltonianLatex"] = SymEngine::latex(*hmexp);
    std::string hs = hmString.substr(0, hmString.size() - 3);
    hs = replaceAll(hs, "**", "^");
    hs.erase(std::remove(hs.begin(), hs.end(), '*'), hs.end());
    result["hamiltonian"] = hs;

    result["portHamiltonianMatrices"] = phs;
    result["success"] = true;
  } catch (BGException &ex) {
    result["error"] = ex.what();
  } catch (std::runtime_error &e) {
    result["error"] = e.what();
  }
  // std::cout << result.dump(2) << std::endl;
  return result;
}