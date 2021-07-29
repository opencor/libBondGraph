std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>>> BondGraph::computeStateEquation()
{
    SymEngine::vec_basic dofs;
    // dx,e,f,x,u
    SymEngine::vec_basic dx;
    SymEngine::vec_basic e;
    SymEngine::vec_basic f;
    SymEngine::vec_basic x;
    SymEngine::vec_basic u;
    std::unordered_map<std::string, std::tuple<std::string,std::string,char> > dimensions; //varible name, value, state, control or parameter
    std::unordered_map<size_t, SymEngine::map_basic_basic> portHamiltonianCoordinates;
    // Global dof matrix
    //Count the number of element in each group
    mUcount = 0; // Number of source bonds
    mScount = 0; // Number of storage bonds with integral causality
    mRcount = 0; // Number of dissipative bonds
    mJcount = 0; // Number of junction bonds
    mSources.clear();
    //Find all the connected components

    std::vector<RCPLIB::RCP<BondGraphElementBase>> connectedComponents;
    std::map<std::string, RCPLIB::RCP<BondGraphElementBase>> compIDMap;
    for (auto &bd : mBonds) {
        auto from = bd->getFromPort()->getComponent();
        auto to = bd->getToPort()->getComponent();
        compIDMap[from->getId()] = from;
        compIDMap[to->getId()] = to;
    }
    for (auto &c : mComponents) { //Maintain order
        auto id = c->getId();
        if (compIDMap.find(id) != compIDMap.end()) {
            //std::cerr << compIDMap[id]->getId() << "\t" << compIDMap[id]->getName() << std::endl;
            connectedComponents.push_back(compIDMap[id]);
        }
    }

    logDebug(" Connected components ", connectedComponents.size(), " num bonds ", mBonds.size());
    long int dofID = 0;
    int portID = 0;
    std::ostringstream ss, dss, ess, fss;

    auto getPreciseUnitString = [](units::precise_unit &un) {
        if (un != units::precise::one) {
            auto mult = un.multiplier();
            auto baseU = un.base_units();
            return units::to_string(units::precise_unit(baseU, mult));
        }
        return units::to_string(un);
    };

    auto getPreciseUnitStringPerSecond = [](units::precise_unit &un) {
        if (un != units::precise::one) {
            auto mult = un.multiplier();
            auto baseU = un.base_units();
            return units::to_string(units::precise_unit(baseU, mult) / units::precise::second);
        }
        return units::to_string(units::precise::hertz);
    };

    // Do storage, dissipative, junction, source
    for (auto &mc : connectedComponents) {
        if (mc->getComponentGroup() == eS) {
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
            ess << "e_" << portID;
            fss << "f_" << portID;
            e.push_back(SymEngine::symbol(ess.str()));
            f.push_back(SymEngine::symbol(fss.str()));
            mc->getPorts(0)->setDofIndex(portID);

            auto un = std::get<1>(values[0])->units;
            auto vn = std::get<1>(values[0])->value;
            //auto mult = un.multiplier();
            //auto baseU = un.base_units();
            dimensions[ss.str()] = std::make_tuple(getPreciseUnitString(un),vn,'s'); //units::to_string(units::precise_unit(baseU, mult));
            //logInfo("Dimensions for state ", ss.str(), " = ", std::get<0>(dimensions[ss.str()]));
            dimensions[dss.str()] = std::make_tuple(getPreciseUnitStringPerSecond(un),"0",'d'); //units::to_string(units::precise_unit(baseU, mult) / units::precise::second);

            ++portID;
            ++dofID;
            ++mScount;

            //Dimensions for the parameters
            for (int i = mc->getNumStates(); i < values.size(); i++) {
                //Use the actual name instead of prefix
                std::string pname = std::get<1>(values[i])->name;
                if (dimensions.find(pname) == dimensions.end()) {
                    auto un = std::get<1>(values[i])->units;
                    auto vn = std::get<1>(values[i])->value;
                    dimensions[pname] = std::make_tuple(getPreciseUnitString(un),vn,'p');
                    //logInfo("Dimensions for parameter ", pname, " = ", std::get<0>(dimensions[pname]));
                }
            }
        }
        if (mc->getComponentGroup() == ePH) {
            auto values = mc->values();
            SymEngine::map_basic_basic globalCoordinates;
            //The entire ph has one dofID
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
                ess << "e_" << portID;
                fss << "f_" << portID;
                e.push_back(SymEngine::symbol(ess.str()));
                f.push_back(SymEngine::symbol(fss.str()));
                mc->getPorts(st)->setDofIndex(portID);

                auto un = std::get<1>(values[st])->units;
                auto vn = std::get<1>(values[st])->value;

                dimensions[ss.str()] = std::make_tuple(getPreciseUnitString(un), vn, 's');
                //logInfo("Dimensions for state ", ss.str(), " = ", std::get<0>(dimensions[ss.str()]));
                dimensions[dss.str()] = std::make_tuple(getPreciseUnitStringPerSecond(un), "0", 'd');
                globalCoordinates[SymEngine::parse(stateName)] = SymEngine::parse(ss.str());
                globalCoordinates[SymEngine::parse("dot_" + stateName)] = SymEngine::parse("dot_" + ss.str());
                globalCoordinates[SymEngine::parse("e_" + std::to_string(st))] = SymEngine::parse(ess.str());
                globalCoordinates[SymEngine::parse("f_" + std::to_string(st))] = SymEngine::parse(fss.str());

                ++portID;
                ++mScount;
            }
            //Dimensions for the parameters
            for (int i = mc->getNumStates(); i < values.size(); i++) {
                //Use the actual name instead of prefix
                std::string pname = std::get<1>(values[i])->name;
                globalCoordinates[SymEngine::parse(std::get<0>(values[i]))] = SymEngine::parse(pname);
                if (dimensions.find(pname) == dimensions.end()) {
                    auto un = std::get<1>(values[i])->units;
                    auto vn = std::get<1>(values[i])->value;
                    dimensions[pname] = std::make_tuple(getPreciseUnitString(un), vn, 'p');
                    //logInfo("Dimensions for parameter ", pname, " = ", std::get<0>(dimensions[pname]));
                }
            }
            portHamiltonianCoordinates[dofID] = globalCoordinates;
        }
    }
    for (auto &mc : connectedComponents) {
        if (mc->getComponentGroup() == eR) {
            auto ports = mc->getPorts();
            //Handle reactions
            for (int pi = 0; pi < ports.size(); pi++) {
                ess.str("");
                ess.clear();
                fss.str("");
                fss.clear();
                ess << "e_" << portID;
                fss << "f_" << portID;
                e.push_back(SymEngine::symbol(ess.str()));
                f.push_back(SymEngine::symbol(fss.str()));
                ports[pi]->setDofIndex(portID);
                ++portID;
            }
            mc->setDof(dofID);
            ++dofID;
            ++mRcount;

            auto values = mc->values();
            for (int i = mc->getNumStates(); i < values.size(); i++) {
                //Use the actual name instead of prefix
                std::string pname = std::get<1>(values[i])->name;
                if (dimensions.find(pname) == dimensions.end()) {
                    auto un = std::get<1>(values[i])->units;
                    auto vn = std::get<1>(values[i])->value;
                    dimensions[pname] = std::make_tuple(getPreciseUnitString(un), vn, 'p');
                    //logInfo("Dimensions for parameter ", pname, " = ", std::get<0>(dimensions[pname]));
                }
            }
        }
    }
    for (auto &mc : connectedComponents) {
        if (mc->getComponentGroup() == eU) {
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
            auto mcValues = mc->values();
            for (int pi = mc->getNumStates(); pi < mcValues.size(); pi++) {
                ess.str("");
                ess.clear();
                ess << std::get<0>(mcValues[pi]) << "_" << dofID;
                u.push_back(SymEngine::symbol(ess.str()));

                auto un = std::get<1>(mcValues[pi])->units;
                auto vn = std::get<1>(mcValues[pi])->value;

                dimensions[ess.str()] =  std::make_tuple(getPreciseUnitString(un),vn,'c');
                //logInfo("Dimensions for state ", ess.str(), " = ", std::get<0>(dimensions[ess.str()]));
                //Control variables have a derivative term
                dimensions["dot_" + ess.str()] = std::make_tuple(getPreciseUnitStringPerSecond(un),"0",'c');
            }
            mSources.push_back(mc);
            ++portID;
            ++dofID;
            ++mUcount;
        }
    }

    for (auto &mc : connectedComponents) {
        if (mc->getComponentGroup() == eJ) {
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
                ports[i]->setDofIndex(portID);
                ++portID;
            }
            mc->setDof(dofID);
            ++dofID;
            ++mJcount;

            if (mc->getType() != eZero && mc->getType() != eOne) { //Transformers are junctions too, so check specifically for 0-,1- junctions
                auto values = mc->values();
                for (int i = mc->getNumStates(); i < values.size(); i++) {
                    //Use the actual name instead of prefix
                    std::string pname = std::get<1>(values[i])->name;
                    if (dimensions.find(pname) == dimensions.end()) {
                        auto un = std::get<1>(values[i])->units;
                        auto vn = std::get<1>(values[i])->value;
                        dimensions[pname] = std::make_tuple(getPreciseUnitString(un), vn, 'p');
                        //logInfo("Dimensions for parameter ", pname, " = ", std::get<0>(dimensions[pname]));
                    }
                }
            }
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

    //JS matrix
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
    //Handle constitutive relations
    //Create coordinate map
    coordinateMap.clear();
    long int di = 0;
    for (auto &c : dofs) {
        ss.str("");
        ss.clear();
        ss << *c;
        coordinateMap[ss.str()] = di++;
    }
    std::unordered_map<unsigned int, std::vector<ExpressionTermsMap>> cIndexes;
    for (auto &mc : connectedComponents) {
        std::vector<ExpressionTermsMap> cx = getLinearAndNonlinearTerms(mc, coordinateMap);
        cIndexes[rows] = cx;
        for (int rc = 0; rc < mc->getConstitutiveEquations().size(); rc++) {
            mc->setConstitutiveEqIndex(rows, rc);
            rows++;
        }
    }

    //logDebug("Coordinates");
    //logDebug(coordinates);
    SymEngine::DenseMatrix linOp(rows, dofs.size());
    SymEngine::DenseMatrix nonlinearTerms(rows, 1);
    zeros(nonlinearTerms); //Initialize the matices with zeros
    zeros(linOp);
    for (int i = 0; i < ix.size(); i++) {
        linOp.set(ix[i], iy[i], SymEngine::integer(dir[i]));
    }
    //Constitutive equations
    for (auto &mp : cIndexes) {
        auto rowIx = mp.first;
        for (auto &mx : mp.second) {
            for (int i = 0; i < mx.indexes.size(); i++) {
                auto col = mx.indexes[i];
                auto cef = mx.coefficients[i];
                if (col > -1) {
                    linOp.set(rowIx, col, cef);
                } else {
                    nonlinearTerms.set(rowIx, 0, SymEngine::expand(SymEngine::add(nonlinearTerms.get(rowIx, 0), cef)));
                }
            }
            rowIx++;
        }
    }

    /*
    logDebug("Coordinates ");
    logDebug(coordinates);
    logDebug("Linear terms ");
    logDebug(linOp);
    logDebug(" NonLinear Op ");
    logDebug(nonlinearTerms);
	*/
    auto snf = getSmithNormalForm(linOp, nonlinearTerms);
    //Add constraints from Port Hamiltionians
    for (auto &mc : connectedComponents) {
        if (mc->getComponentGroup() == ePH) {
            RCPLIB::RCP<PortHamiltonian> ph = RCPLIB::rcp_dynamic_cast<PortHamiltonian>(mc);
            const std::vector<std::string> &cons = ph->getConstraints();
            for (auto &c : cons) {
                //Go from local to global coordinates
                auto cexp = SymEngine::parse(c);
                long int id = ph->getDof();
                auto &subs = portHamiltonianCoordinates[id];
                snf.constraints.push_back(cexp->subs(subs));
            }
        }
    }
    /*   
    logDebug(" Linear Op ");
    logDebug(snf.linearOp);
    logDebug(" NonLinear Op ");
    logDebug(snf.nonlinearOp);
    logDebug(" Constraints ");
    logDebug(snf.constraints);
	*/

    auto subsExprs = findSubstitutions(snf.linearOp, snf.nonlinearOp, snf.constraints, coordinates, numStates, numBonds);
    substituteValues(snf.nonlinearOp, subsExprs);
    substituteValues(snf.constraints, subsExprs);

    //Process constraints
    auto sn = process_constraints(snf, coordinates, numStates, numBonds, numControlVariables);
    auto subsExprsC = findSubstitutions(sn.linearOp, sn.nonlinearOp, sn.constraints, coordinates, numStates, numBonds);
    substituteValues(sn.nonlinearOp, subsExprsC);
    substituteValues(sn.constraints, subsExprsC);

    SymEngine::DenseMatrix &linearOp = sn.linearOp;
    SymEngine::DenseMatrix &nonlinearOp = sn.nonlinearOp;
    SymEngine::vec_basic &constraints = sn.constraints;

    int rows_added = 0;
    std::vector<int> added_cvs;
    SymEngine::map_basic_basic cv_diff_dict;
    std::vector<std::tuple<int, int, SymEngine::RCP<const SymEngine::Basic>>> lin_dict;

    //Create coordinate map
    coordinateMap.clear();
    di = 0;
    for (auto &c : dofs) {
        ss.str("");
        ss.clear();
        ss << *c;
        std::string dofname = ss.str(); //Eliminate memory loss issue
        coordinateMap[dofname] = di++;
    }
    bool solvable = true; //Flag to check if the equations are solvable.. Some incorrect bg formulations will lead to state and bond variables being solved to zero
    // - Linear constraints; ie Lx = 0
    // - Nonlinear Constraints Lx + F(x) = 0
    //
    // Linear constraints are rows with more than 1 non-zero
    // that are not in the derivative subspace, and have a zero nonlinear part
    offset = 2 * numBonds + numStates;
    if (linearOp.nrows() < offset) {
        SymEngine::set_basic c_atoms;
        for (int i = 0; i < coordinates.nrows(); i++) {
            auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*coordinates.get(i, 0));
            c_atoms.insert(atoms.begin(), atoms.end());
        }
        for (int row = offset - 1; row > linearOp.nrows(); row--) {
            auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*nonlinearOp.get(row, 0));
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

    std::vector<std::tuple<int, int, SymEngine::RCP<const SymEngine::Basic>>> cv_dict;
    rows = linearOp.nrows(); //Only parse for existing rows, not the ones that will be added
    for (int row = offset; row < rows; row++) {
        auto nonlinear_constraint = nonlinearOp.get(row, 0);
        //Size change during the call so do not allocate outside
        SymEngine::DenseMatrix linearOpSub(1, linearOp.ncols() - offset);
        linearOp.submatrix(linearOpSub, row, offset, row, linearOp.ncols() - 1);
        if (is_true(linearOpSub.is_zero()) && SymEngine::eq(*nonlinear_constraint, *SymEngine::zero)) {
            continue;
        }
        //Size change during the call so do not allocate outside
        SymEngine::DenseMatrix state_constraint(1, numStates);
        SymEngine::DenseMatrix control_constraint(1, linearOp.ncols() - offset - numStates);
        SymEngine::DenseMatrix sj(1, offset + numControlVariables);

        linearOp.submatrix(state_constraint, row, offset, row, offset + numStates - 1);
        linearOp.submatrix(control_constraint, row, offset + numStates, row, linearOp.ncols() - 1);
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
                if (std::find(added_cvs.begin(), added_cvs.end(), cv_col) == added_cvs.end()) {
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
                    idx = std::find(added_cvs.begin(), added_cvs.end(), cv_col) - added_cvs.begin();
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
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(ss.str());
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
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(ss.str());
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
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(ss.str());
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
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(ss.str());
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
                sumxy = SymEngine::add(sumxy, SymEngine::mul(jac_x[x], coordinates.get(x, 0)));
            }
            auto expr = SymEngine::expand(sumxy);
            SymEngine::RCP<const SymEngine::Basic> num, denom;
            SymEngine::as_numer_denom(expr, SymEngine::outArg(num), SymEngine::outArg(denom));
            if (is_true(state_constraint.is_zero())) {
                std::map<long int, SymEngine::RCP<const SymEngine::Basic>> ld;
                SymEngine::RCP<const SymEngine::Basic> nl;
                SymEngine::map_basic_basic dummy;
                ss.str("");
                ss.clear();
                ss << *num;
                std::string numer = ss.str();
                std::tie(ld, nl) = getLinearCoefficientsAndNonlinearTerms(numer, coordinateMap, dummy);
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

    auto System = SymEngine::DenseMatrix(linearOp.nrows(), coordinates.ncols());
    mul_dense_dense(linearOp, coordinates, System);

    SymEngine::map_basic_basic equations;
    SymEngine::map_basic_basic bondEquations;
    SymEngine::map_basic_basic bondSubs;

    /*
    //Load what is know by definition in the components
    for(auto&  elem : connectedComponents) {
        std::cerr<<__LINE__<<"\t Element name "<<elem->getName()<<std::endl;
        if (elem->getComponentGroup() == eU) {
            //Has control variable(s)
            auto values = elem->values();
            for (int i = elem->getNumStates(); i < values.size(); i++) {
                ss.str("");
                ss.clear();
                std::string controlVar = std::get<0>(values[i]);
                ss << controlVar << "_" << elem->getId();
                auto un = std::get<1>(values[i])->units;
                auto mult = un.multiplier();
                auto baseU = un.base_units();
                dimensions[controlVar] = units::to_string(units::precise_unit(baseU, mult));
            }
        }
        if (elem->getComponentGroup() == eS) {
            //Has a state variable(s)
            auto values = elem->values();
            for (int i = 0; i < elem->getNumStates(); i++) {
                ss.str("");
                ss.clear();

                std::string stateName = std::get<0>(values[i]);
                auto ploc = stateName.rfind("_");
                std::string statePrefix = stateName;
                if (ploc != std::string::npos) {
                    statePrefix = stateName.substr(0, ploc + 1);
                    ss << statePrefix << elem->getId();
                } else {
                    ss << stateName << elem->getId();
                }
                std::string stateVar = ss.str();
                std::string dstateVar = "dot_" + stateVar;
                //ss.str("");
                //ss.clear();
                //ss << std::get<0>(values[i]) << elem->getId();
                auto un = std::get<1>(values[i])->units;
                auto mult = un.multiplier();
                auto baseU = un.base_units();
                dimensions[stateVar] = units::to_string(units::precise_unit(baseU, mult));
                dimensions[dstateVar] = units::to_string(units::precise_unit(baseU, mult) / units::precise::second);
                std::cerr<<__LINE__<<"\t"<<stateVar<<"\t"<<dimensions[stateVar]<<std::endl;
            }
        }
        if (elem->getType() != eZero && elem->getType() != eOne) { //Transformers are junctions too, so check specifically for 0-,1- junctions
            auto values = elem->values();
            for (int i = elem->getNumStates(); i < values.size(); i++) {
                ss.str("");
                ss.clear();
                //Use the actual name instead of prefix
                std::string pname = std::get<1>(values[i])->name;
                if (dimensions.find(pname) == dimensions.end()) {
                    auto un = std::get<1>(values[i])->units;
                    auto mult = un.multiplier();
                    auto baseU = un.base_units();
                    dimensions[pname] = units::to_string(units::precise_unit(baseU, mult));
                    std::cerr<<__LINE__<<"\t"<<pname<<"\t"<<dimensions[pname]<<std::endl;
                }
            }
        }
    }*/

    for (int i = numStates; i < numStates + 2 * numBonds; i++) {
        ss.str("");
        ss.clear();
        auto lhs = System.get(i, 0);
        auto rhs = nonlinearOp.get(i, 0);
        if (!SymEngine::eq(*lhs, *SymEngine::zero) || !SymEngine::eq(*rhs, *SymEngine::zero)) {
            ss << *lhs << " + " << *rhs;
            bondSubs[lhs] = rhs;
            auto expr = SymEngine::parse(ss.str());
            //Expression will have both the coordinate and related terms, as lhs is not always the coordinate values
            //Solve and get the result for the coordinate
            ss.str("");
            ss.clear();
            ss << *coordinates.get(i, 0);
            auto variable = ss.str();
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(variable);
            auto solns = SymEngine::solve(expr, csym)->get_args();
            if (SymEngine::eq(*solns[0], *SymEngine::zero)) {
                logWarn(*lhs, " + ", *rhs, " gives solution of 0 for ", variable);
                solvable = false;
            }
            bondEquations[coordinates.get(i, 0)] = solns[0];

            //Find the dimension of the varible if not already found
            if (dimensions.find(variable) == dimensions.end()) {
                //Parse through the terms in the solution and determine the dimensions of each term
                auto ed = getDimensions(solns[0], dimensions, 'v');
                //std::cerr << __LINE__ << "\t" << *expr << "\t" << ed << " for " << *solns[0] << std::endl;
                //logInfo("Dimensions for ", variable, " = ", *solns[0], " \t ", std::get<0>(ed));
                dimensions[variable] = ed;
            }
        }
    }

    for (int i = 0; i < numStates; i++) {
        ss.str("");
        ss.clear();
        auto lhs = System.get(i, 0);
        auto rhs = nonlinearOp.get(i, 0);
        if (!SymEngine::eq(*lhs, *SymEngine::zero) || !SymEngine::eq(*rhs, *SymEngine::zero)) {
            ss << *lhs << " + " << *rhs;
            auto expr = SymEngine::parse(ss.str());
            //Expression will have both the coordinate and related terms, as lhs is not always the coordinate values
            //Solve and get the result for the coordinate
            ss.str("");
            ss.clear();
            ss << *coordinates.get(i, 0);
            auto variable = ss.str();
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(variable);
            auto solns = SymEngine::solve(expr, csym)->get_args();
            if (SymEngine::eq(*solns[0], *SymEngine::zero)) {
                logWarn(*lhs, " + ", *rhs, " gives solution of 0 for ", variable);
                solvable = false;
            }
            equations[coordinates.get(i, 0)] = solns[0];

            //Find the dimension of the varible if not already found
            if (dimensions.find(variable) == dimensions.end()) {
                //Parse through the terms in the solution and determine the dimensions of each term
                auto ed = getDimensions(solns[0], dimensions,'d');
                //logInfo("Dimensions for ", variable, " = ", *solns[0], " \t ", std::get<0>(ed));
                dimensions[variable] = ed;
            }
        }
    }

    for (int i = numStates + 2 * numBonds; i < System.nrows(); i++) {
        ss.str("");
        ss.clear();
        auto lhs = System.get(i, 0);
        auto rhs = nonlinearOp.get(i, 0);
        if (!SymEngine::eq(*lhs, *SymEngine::zero) || !SymEngine::eq(*rhs, *SymEngine::zero)) {
            ss << *lhs << " + " << *rhs;
            auto expr = SymEngine::parse(ss.str());
            //Expression will have both the coordinate and related terms, as lhs is not always the coordinate values
            //Solve and get the result for the coordinate
            ss.str("");
            ss.clear();
            ss << *coordinates.get(i, 0);
            auto variable = ss.str();
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(variable);
            auto solns = SymEngine::solve(expr, csym)->get_args();
            if (SymEngine::eq(*solns[0], *SymEngine::zero)) {
                //std::cout << *lhs << " + " << *rhs << " gives solution of 0 for " << variable << std::endl;
                logWarn(*lhs, " + ", *rhs, " gives solution of 0 for ", variable);
                solvable = false;
            }
            equations[coordinates.get(i, 0)] = solns[0];

            //Find the dimension of the varible if not already found
            if (dimensions.find(variable) == dimensions.end()) {
                //Parse through the terms in the solution and determine the dimensions of each term
                auto ed = getDimensions(solns[0], dimensions,'v');
                //logInfo("Dimensions for ", variable, " = ", *solns[0], " \t ", std::get<0>(ed));
                dimensions[variable] = ed;
            }
        }
    }
    //parse constraints for state derivatives
    SymEngine::vec_basic minConstraints;
    for (int i = 0; i < constraints.size(); i++) {
        auto cons = constraints[i];
        bool removeCons = false;
        for (int j = 0; j < numStates; j++) {
            ss.str("");
            ss.clear();
            ss << *coordinates.get(j, 0);
            SymEngine::RCP<const SymEngine::Symbol> csym = SymEngine::symbol(ss.str());
            auto solns = SymEngine::solve(cons, csym)->get_args();
            if (solns.size() > 0) {
                equations[coordinates.get(j, 0)] = solns[0];
                auto variable = ss.str();
                //Find the dimension of the varible if not already found
                if (dimensions.find(variable) == dimensions.end()) {
                    //Parse through the terms in the solution and determine the dimensions of each term
                    auto ed = getDimensions(solns[0], dimensions,'v');
                    //logInfo("Dimensions for ", variable," = ", *solns[0] ," \t ", std::get<0>(ed));
                    dimensions[variable] = ed;
                }
                removeCons = true;
                break;
            }
        }
        if (!removeCons)
            minConstraints.push_back(cons);
    }

    //std::cerr << equations << std::endl;
    computedCoordinates = coordinates;
    //std::cerr<<bondEquations<<std::endl;
    return std::make_tuple(solvable, equations, bondEquations, minConstraints, dimensions);
}