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
#include "Bond.h"
#include "Elements.h"
#include "Elementsbase.h"
#include "Port.h"
#include <algorithm>
#include <exception>
#include <set>
#include <sstream>
#include <symengine/eval.h>
#include <symengine/matrix.h>
#include <symengine/parser.h>
#include <symengine/parser/parser.h>
#include <symengine/solve.h>
#include <symengine/visitor.h>
#include <symengine/simplify.h>
#include <tuple>
#include <units.hpp>
#include <unordered_map>
#include <vector>
//#include <regex>
//#include <boost/regex.hpp>
#include "simplifyexplog.h"


/**
 * @brief Collection of algebraic manipulation routines used to determine the reduced system of equations
 * 
 * Not abstracted into a class for ease of unit testing
 */

namespace BG {

/**
     * @brief Struture to hold linear and nonlinear terms
     * 
     * indexes member has a non zero index into the linear operator, when this index is -1 the corresponding term in coefficients is a nonlinear term
     * coefficients member stores the symbolic coefficient or nonlinear term corresponding to the indexes entry
     */
class ExpressionTermsMap
{
public:
    std::vector<long int> indexes;
    SymEngine::vec_basic coefficients;

    ExpressionTermsMap(std::vector<long int> &ix, SymEngine::vec_basic &c)
        : indexes(ix)
        , coefficients(c)
    {
    }
};

/**
     * @brief Structure to encapsulate the results from Smith Normal Form calculation
     * linearOp - the linear operator in upper triangular form with diagonal entries being 1
     * nonlinearOp - vector of nonlinear terms
     * constraints - vector of constraints
     */
class SmithNormalForm
{
public:
    SymEngine::DenseMatrix linearOp;
    SymEngine::DenseMatrix nonlinearOp;
    SymEngine::vec_basic constraints;

    SmithNormalForm(SymEngine::DenseMatrix &l, SymEngine::DenseMatrix &n, SymEngine::vec_basic &c)
        : linearOp(l)
        , nonlinearOp(n)
        , constraints(c)
    {
    }
};

/**
     * @brief Get the coefficients for the linear terms and the nonlinear terms associated with a constitutive equations
     * 
     * @param ceq - Constitutive equation (std:string)
     * @param coordinates - unordered map from state,port variable name, or other symbolic name with entry in the global dof to index into that global dof array
     * @param usedNames - map of symbols used by the element and their global naming, includes parameters and names of port variables associated with the bondgraph element to which the constitutive equation belongs
     * @return std::tuple<std::map<long int,SymEngine::RCP<const SymEngine::Basic>>, SymEngine::RCP<const SymEngine::Basic> > 
     */
std::tuple<std::map<long int, SymEngine::RCP<const SymEngine::Basic>>, SymEngine::RCP<const SymEngine::Basic>> getLinearCoefficientsAndNonlinearTerms(std::string &ceq,
                                                                                                                                                      std::unordered_map<std::string, long int> &coordinates,
                                                                                                                                                      SymEngine::map_basic_basic &usedNames)
{
    std::map<long int, SymEngine::RCP<const SymEngine::Basic>> coeff_dict; //Dictionary to store linear terms, dof index -> term/coefficient
    SymEngine::RCP<const SymEngine::Basic> nonlinear_terms = SymEngine::parse("0");
    auto equation = SymEngine::parse(ceq);
    /*
            Split the equation into additive terms
            parse each additive term to find coefficients of the term if the term multiple of a dof entry
            enter a item in coeff_dict for the dof and coefficient
            if the term has a dof entry in the denominator or as a argument to a function or is multiplied 
            with another dof entry, then report this term as a nonlinear term
            In case there are more than one nonlinear term, they are summed.
        */
    auto terms = SymEngine::expand(equation)->get_args();
    std::ostringstream ss;
    for (auto &term : terms) {
        SymEngine::vec_basic factors; // expand term as arguments of a product, in case of an non product term multiply this by 1
        if (SymEngine::is_a<SymEngine::Mul>(*term)) {
            auto comp = term->get_args();
            factors.insert(factors.end(), comp.begin(), comp.end());
        } else {
            factors.push_back(SymEngine::integer(1));
            factors.push_back(term);
        }
        //Parse for base symbolic term and form its coefficients
        SymEngine::RCP<const SymEngine::Basic> coeff = SymEngine::integer(1);
        SymEngine::vec_basic base;
        for (int fi = factors.size() - 1; fi > -1; fi--) {
            auto c = factors[fi];

            if (SymEngine::is_a<SymEngine::Integer>(*c) || SymEngine::is_a<SymEngine::Rational>(*c)) {
                coeff = SymEngine::mul(c, coeff);
            } else if (SymEngine::is_a<SymEngine::Symbol>(*c) && usedNames.find(c) == usedNames.end()) {
                coeff = SymEngine::mul(c, coeff);
            } else {
                base.push_back(c);
            }
        }

        if (base.size() == 1 && usedNames.find(base[0]) != usedNames.end()) {
            ss.str("");
            ss.clear();
            ss << *usedNames[base[0]];
            coeff_dict[coordinates[ss.str()]] = coeff;
        } else {
            //Change local names to their global names
            SymEngine::RCP<const SymEngine::Basic> new_term = term->subs(usedNames);
            nonlinear_terms = SymEngine::add(new_term, nonlinear_terms);
        }
    }
    return std::make_tuple(coeff_dict, nonlinear_terms);
}

/**
     * @brief Get the Linear And Nonlinear Terms associated with a Bondgraph element (constitutive equations)
     * 
     * @param elem - Bondgraph element of interest
     * @param coordinates - mapping between local state, port element names, and control variables (source variables) to global dof index
     * @return std::vector<ExpressionTermsMap> 
     */
std::vector<ExpressionTermsMap> getLinearAndNonlinearTerms(const RCPLIB::RCP<BGElement>  &elem_, std::unordered_map<std::string, long int> &coordinates)
{
    const RCPLIB::RCP<BondGraphElementBase>  elem = RCPLIB::rcp_dynamic_cast<BondGraphElementBase>(elem_);
    std::vector<ExpressionTermsMap> mapping;
    std::vector<std::string> &lceq = elem->getConstitutiveEquations();
    std::unordered_map<std::string, std::string> usedNames; //Local to global dof name map
    SymEngine::map_basic_basic usedNamesSymbol;
    std::unordered_map<std::string, std::string> parameters;
    //Update parameter names in constitutive equations to match global names
    auto values = elem->values();
    std::vector<std::string> ceq(lceq);
    for (int ix = elem->getNumStates(); ix < values.size(); ix++) {
        auto target = std::get<1>(values[ix]);
        if (target->universalConstant)
            continue;
        for (int i = 0; i < ceq.size(); i++) {
            ceq[i] = replaceAll(ceq[i], target->prefix, target->name);
        }
    }

    //Create the element specific usednames map, the usedname's value will have an entry in coordinates
    std::ostringstream ss;
    //Handle sources
    if (elem->getComponentGroup() == eU) {
        //Has control variable(s)
        for (int i = elem->getNumStates(); i < values.size(); i++) {
            ss.str("");
            ss.clear();
            std::string controlVar = std::get<0>(values[i]);
            ss << controlVar << "_" << elem->getDof();
            std::string dofName = ss.str();
            //usedNames[controlVar] = dofName; //Incorrect
            //Control variable names in construtive equations are updated in section [Update parameter names]
            usedNames[dofName] = dofName;
        }
    }
    if (elem->getComponentGroup() == eS) {
        //Has a state variable(s)
        for (int i = 0; i < elem->getNumStates(); i++) {
            ss.str("");
            ss.clear();
            std::string stateName = std::get<0>(values[i]);
            auto ploc = stateName.rfind("_");
            std::string statePrefix = stateName;
            if (ploc != std::string::npos) {
                statePrefix = stateName.substr(0, ploc + 1);
                ss << statePrefix << i;
            } else {
                ss << stateName << i;
            }
            std::string stateVar = ss.str();
            std::string dstateVar = "dot_" + stateVar;
            ss.str("");
            ss.clear();
            ss << statePrefix << elem->getDof();
            std::string dofName = ss.str();
            usedNames[stateVar] = dofName;
            usedNames[dstateVar] = "dot_" + dofName;
        }
    }
    if (elem->getType() != eZero && elem->getType() != eOne) {
        for (int i = elem->getNumStates(); i < values.size(); i++) {
            ss.str("");
            ss.clear();
            std::string pname = std::get<0>(values[i]);
            if (usedNames.find(pname) == usedNames.end()) {
                parameters[pname] = pname;
            }
        }
    }

    auto ports = elem->getPorts();

    int pi = 0;
    for (auto &p : ports) {
        ss.str("");
        ss.clear();
        ss << "_" << pi;
        std::string ep = "e" + ss.str();
        std::string fp = "f" + ss.str();
        ss.str("");
        ss.clear();
        ss << "_" << p->dofIndex();
        usedNames[ep] = "e" + ss.str();
        usedNames[fp] = "f" + ss.str();
        pi++;
    }
    //Create a symbolic version for getLinearCoefficientsAndNonlinearTerms call

    for (auto &un : usedNames) {
        usedNamesSymbol[SymEngine::parse(un.first)] = SymEngine::parse(un.second);
    }

    //Find linear terms coefficients
    //Update nonlinear terms with stateVar names
    std::map<long int, SymEngine::RCP<const SymEngine::Basic>> coeff_dict;
    SymEngine::RCP<const SymEngine::Basic> nonlinear_terms;

    for (auto &eq : ceq) {
        auto expr = SymEngine::parse(eq);
        std::tie(coeff_dict, nonlinear_terms) = getLinearCoefficientsAndNonlinearTerms(eq, coordinates, usedNamesSymbol);
        std::vector<long int> linearIndex;
        SymEngine::vec_basic coeff;
        for (auto &c : coeff_dict) {
            linearIndex.push_back(c.first);
            coeff.push_back(c.second);
        }
        if (!SymEngine::eq(*nonlinear_terms, *SymEngine::zero)) {
            linearIndex.push_back(-1);
            coeff.push_back(nonlinear_terms);
        }
        mapping.push_back(ExpressionTermsMap(linearIndex, coeff));
    }
    return mapping;
}

/**
    * @brief Computes the Smith normal form of the given matrix.
    * 
    * @return SmithNormalForm 
    */
SmithNormalForm getSmithNormalForm(const SymEngine::DenseMatrix &linear, const SymEngine::DenseMatrix &nonlinear)
{
    SymEngine::DenseMatrix M(linear);

    M.row_join(nonlinear);

    auto k = nonlinear.ncols();
    auto r = M.nrows();
    auto c = M.ncols();
    //Compute the reduced row echelon form of M
    SymEngine::permutelist pl;

    SymEngine::vec_basic constraints;
    int pivot = 0;
    auto m = M.ncols() - k;
    for (int col = 0; col < m; col++) {
        if (SymEngine::eq(*M.get(pivot, col), *SymEngine::zero)) {
            int j = -1;
            auto v_max = M.get(pivot, col);
            for (int row = pivot; row < M.nrows(); row++) {
                auto val = M.get(row, col);
                auto v = SymEngine::abs(val);
                if (SymEngine::is_a<SymEngine::Symbol>(*v)) {
                    j = row;
                    v_max = v;
                } else {
                    if (SymEngine::unified_compare(v, v_max) == 1) //! \return -1, 0, 1 for a < b, a == b, a > b
                    {
                        j = row;
                        v_max = v;
                    }
                }
            }
            if (j == -1)
                continue; // all zeros below, skip on to next column
            else
                row_exchange_dense(M, pivot, j);
        }
        auto a = M.get(pivot, col);

        for (int i = 0; i < M.nrows(); i++) {
            if (i != pivot && !SymEngine::eq(*M.get(i, col), *SymEngine::zero)) {
                auto b = SymEngine::div(M.get(i, col), a);
                for (int cc = 0; cc < M.ncols(); cc++) {
                    auto cv = M.get(i, cc);
                    //M.set(i, cc, SymEngine::simplify(SymEngine::sub(cv, SymEngine::mul(b, M.get(pivot, cc)))));
                    M.set(i, cc, SymEngine::sub(cv, SymEngine::mul(b, M.get(pivot, cc))));
                }
            }
        }
        for (int cc = 0; cc < M.ncols(); cc++) {
            auto cv = M.get(pivot, cc);
            //M.set(pivot, cc, SymEngine::simplify(SymEngine::div(cv, a)));
            M.set(pivot, cc, SymEngine::div(cv, a));
        }

        pivot += 1;

        if (pivot >= M.nrows())
            break;
    }

    SymEngine::DenseMatrix Mp = SymEngine::DenseMatrix(c - k, c);
    zeros(Mp); //Populate the matrix else setting or processing will lead to memory faults (its not a SparseMatrix)
    //Temporary matrices
    SymEngine::DenseMatrix subm = SymEngine::DenseMatrix(1, k);
    SymEngine::DenseMatrix suba = SymEngine::DenseMatrix(1, c);

    for (int row = 0; row < r; row++) {
        int leading_coeff = -1;
        for (int col = row; col < c - k; col++) {
            if (!SymEngine::eq(*M.get(row, col), *SymEngine::zero)) {
                leading_coeff = col;
                break;
            }
        }
        //Terms with negative linear coefficient are identified as constraints
        //A constraint is created by summing other dof terms appearing in that row
        if (leading_coeff < 0) {
            M.submatrix(subm, row, c - k, row, c - 1);
            if (!is_true(subm.is_zero())) {
                M.submatrix(suba, row, 0, row, c - 1);
                constraints.push_back(SymEngine::add(suba.as_vec_basic()));
            }
        } else {
            for (int col = 0; col < c; col++) {
                Mp.set(leading_coeff, col, M.get(row, col));
            }
        }
    }
    SymEngine::DenseMatrix linearop = SymEngine::DenseMatrix(c - k, c - k);
    SymEngine::DenseMatrix nonlinearop = SymEngine::DenseMatrix(c - k, k);
    Mp.submatrix(linearop, 0, 0, c - k - 1, c - k - 1);
    Mp.submatrix(nonlinearop, 0, c - k, c - k - 1, c - 1);
    for (int i = 0; i < nonlinearop.nrows(); i++)
        for (int j = 0; j < nonlinearop.ncols(); j++){
            nonlinearop.set(i, j, SymEngine::simplify(nonlinearop.get(i, j)));
        }

    return SmithNormalForm(linearop, nonlinearop, constraints);
}

/**
     * @brief Compute the Intersection between two symbolic sets
     * 
     * @param a 
     * @param b 
     * @return SymEngine::set_basic 
     */
SymEngine::set_basic setIntersection(SymEngine::set_basic &a, SymEngine::set_basic &b)
{
    SymEngine::set_basic res;
    const auto iEnd = a.end();
    for (auto &c : b) {
        if (a.find(c) != iEnd) {
            res.insert(c);
        }
    }
    return res;
}

/**
     * @brief Find solutions for dofs based on current linear operators, nonlinear operators and constraints
     *        Lx + F(x) = 0 is the initial bond graph system
     *        Smith normal form takes this to the form Ix = (I - L)x - F(x) = Rx - F(x)
     *        Here L is the Smith normal form, then if (Rx)_{ii} = 0, and F(x)_i doesn't depend upon x_i
     *        then we have x_i = (Rx)_i - F_i(x)
     * @param linearOp 
     * @param nonlinearOp 
     * @param constraints 
     * @param coordinates 
     * @param numStates 
     * @param numBonds 
     * @return SymEngine::map_basic_basic 
     */
SymEngine::map_basic_basic findSubstitutions(
    SymEngine::DenseMatrix &linearOp,
    SymEngine::DenseMatrix &nonlinearOp,
    SymEngine::vec_basic &constraints,
    SymEngine::DenseMatrix &coordinates,
    unsigned int numStates,
    unsigned int numBonds)
{
    SymEngine::set_basic c_atoms;
    SymEngine::set_basic nonlinear_atoms;
    std::vector<SymEngine::set_basic> nonlinearOpatoms;
    SymEngine::set_basic s_atoms;
    SymEngine::map_basic_basic results;
    for (int i = 0; i < coordinates.nrows(); i++) {
        auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*coordinates.get(i, 0));
        c_atoms.insert(atoms.begin(), atoms.end());
    }
    for (int i = 0; i < nonlinearOp.nrows(); i++) {
        auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*nonlinearOp.get(i, 0));
        nonlinear_atoms.insert(atoms.begin(), atoms.end());
        nonlinearOpatoms.push_back(atoms);
    }
    //Set intersection
    s_atoms = setIntersection(c_atoms, nonlinear_atoms);
    for (auto &cons : constraints) {
        auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*cons);
        auto interS = setIntersection(atoms, c_atoms);
        for (auto &c : interS) {
            s_atoms.insert(c);
        }
    }

    if (s_atoms.size() < 1)
        return results;

    SymEngine::DenseMatrix Rx1(linearOp.nrows(), linearOp.ncols());
    SymEngine::eye(Rx1);
    SymEngine::DenseMatrix Rx(Rx1);
    SymEngine::DenseMatrix negLinearOp(linearOp);
    mul_dense_scalar(linearOp, SymEngine::integer(-1), negLinearOp);
    add_dense_dense(Rx1, negLinearOp, Rx);
    auto loCols = linearOp.ncols();
    SymEngine::DenseMatrix Rxcol(1, loCols);
    SymEngine::DenseMatrix RxmulC(1, 1);
    for (int i = 2 * (numStates + numBonds) - 1; i > -1; i--) {
        auto co = coordinates.get(i, 0);

        if (SymEngine::eq(*Rx.get(i, i), *SymEngine::zero) && s_atoms.find(co) != s_atoms.end() && nonlinearOpatoms[i].find(co) == nonlinearOpatoms[i].end()) {
            Rx.submatrix(Rxcol, i, 0, i, loCols - 1);
            mul_dense_dense(Rxcol, coordinates, RxmulC);
            auto eqn = SymEngine::sub(RxmulC.get(0, 0), nonlinearOp.get(i, 0));
            if (!SymEngine::eq(*eqn, *SymEngine::zero)) {
                SymEngine::map_basic_basic lres;
                lres[co] = eqn;
                for (auto &c : results) {
                    results[c.first] = c.second->subs(lres);
                }
                results[co] = eqn;
            }
        }
    }
    return results;
}

/**
     * @brief Substitute expressions in 'substitutions' dict in to the expression 'expr'
     * 
     * @param expr 
     * @param substitutions 
     * @return SymEngine::RCP<const SymEngine::Basic> 
     */
inline SymEngine::RCP<const SymEngine::Basic> substituteValues(SymEngine::RCP<const SymEngine::Basic> &expr, SymEngine::map_basic_basic &substitutions)
{
    auto sexpr = SymEngine::expand(expr->subs(substitutions));
    std::ostringstream ss;
    //Hack, as symengine currently (May 2021) does not provide simplify
    //Replace log(exp(*)) or exp(log(*)) by ()
    ss << *sexpr;
    auto expString = ss.str();
    auto nExpr = SymEngine::simplify(SymEngine::parse(expString));
    return nExpr;
}

/**
     * @brief Substitute expressions in 'substitutions' dict in to each expression in matrix
     * 
     * @param matrix 
     * @param substitutions 
     */
void substituteValues(SymEngine::DenseMatrix &matrix, SymEngine::map_basic_basic substitutions)
{
    for (int i = 0; i < matrix.nrows(); i++) {
        for (int j = 0; j < matrix.ncols(); j++) {
            auto expr = matrix.get(i, j);
            if (!SymEngine::eq(*expr, *SymEngine::zero)) {
                matrix.set(i, j, substituteValues(expr, substitutions));
            }
        }
    }
}

/**
     * @brief Substitute expressions in 'substitutions' dict in to each expression in vec
     * 
     * @param vec 
     * @param substitutions 
     */
void substituteValues(SymEngine::vec_basic &vec, SymEngine::map_basic_basic substitutions)
{
    for (int i = 0; i < vec.size(); i++) {
        auto expr = vec[i];
        if (!SymEngine::eq(*expr, *SymEngine::zero)) {
            vec[i] = substituteValues(expr, substitutions);
        }
    }
}

/**
    * @brief Find the fractional form of an expression
    * 
    * @param expr 
    * @return SymEngine::vec_basic numerator, denominator
    */
SymEngine::vec_basic fraction(SymEngine::RCP<const SymEngine::Basic> &expr)
{
    SymEngine::RCP<const SymEngine::Basic> num, denom;
    SymEngine::as_numer_denom(expr, SymEngine::outArg(num), SymEngine::outArg(denom));
    SymEngine::vec_basic res;
    res.push_back(num);
    res.push_back(denom);
    return res;
}

/**
     * @brief Process the constraints and reduce their order where possible
     * 
     * @param ops 
     * @param coordinates 
     * @param numStates 
     * @param numBonds 
     * @param numControlVars 
     * @return SmithNormalForm 
     */
SmithNormalForm process_constraints(SmithNormalForm &ops,
                                    SymEngine::DenseMatrix &coordinates,
                                    unsigned int numStates,
                                    unsigned int numBonds,
                                    unsigned int numControlVars)
{
    SymEngine::DenseMatrix &linear_op = ops.linearOp;
    SymEngine::DenseMatrix &nonlinear_op = ops.nonlinearOp;
    SymEngine::vec_basic &constraints = ops.constraints;

    SymEngine::vec_basic initial_constraints;
    auto ss_size = numStates;
    auto js_size = numBonds;
    auto cs_size = numControlVars;
    auto n = coordinates.nrows();
    auto numCoords = n - cs_size;
    auto offset = 2 * js_size + ss_size;
    SymEngine::DenseMatrix coord_atoms(offset + ss_size, 1);
    coordinates.submatrix(coord_atoms, 0, 0, offset + ss_size - 1, 0); 

    SymEngine::DenseMatrix linOpCof(n - (offset + ss_size), linear_op.ncols());
    SymEngine::DenseMatrix nonlinOpCof(n - (offset + ss_size), 1);
    SymEngine::DenseMatrix nonlinOpCoord(n - (offset + ss_size), 1);
    SymEngine::DenseMatrix cv_constraints(n - (offset + ss_size), 1);
    //Process control variables
    //Check that these variables exist, else submatrix will fail
    if (offset + ss_size < n) {
        linear_op.submatrix(linOpCof, offset + ss_size, 0, n - 1, linear_op.ncols() - 1);
        nonlinear_op.submatrix(nonlinOpCof, offset + ss_size, 0, n - 1, 0);
        linOpCof.mul_matrix(coordinates, nonlinOpCoord);
        zeros(cv_constraints); //Allocate here when we know space is required
        nonlinOpCoord.add_matrix(nonlinOpCof, cv_constraints);

        for (int i = 0; i < cv_constraints.nrows(); i++) {
            auto cv = cv_constraints.get(i, 0);
            if (!SymEngine::eq(*cv, *SymEngine::zero)) {
                constraints.push_back(cv);
            }
        }
    }

    SymEngine::DenseMatrix linOp(offset + ss_size, linear_op.ncols());
    SymEngine::DenseMatrix nonlinOp(offset + ss_size, 1);
    linear_op.submatrix(linOp, 0, 0, offset + ss_size - 1, linear_op.ncols() - 1);
    nonlinear_op.submatrix(nonlinOp, 0, 0, offset + ss_size - 1, 0);
    if (constraints.size() > 0) {
        SymEngine::set_basic c_atoms; //Set of all coordinate atoms
        for (int i = 0; i < numCoords; i++) {
            auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*coordinates.get(i, 0));
            c_atoms.insert(atoms.begin(), atoms.end());
        }
        std::ostringstream ss;
        for (auto &cons : constraints) {
            SymEngine::vec_basic frac = fraction(cons);
            auto constraint = frac[0];
            //Find all dofs in the constraint
            auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*constraint);
            auto s_atoms = setIntersection(c_atoms, atoms);
            //If there is just one dof, check if there is a solution
            //if there is augment the linear operator and coordinates vector
            if (s_atoms.size() == 1) {
                //Dont know how to convert basic to symbol
                auto ab = *(s_atoms.begin());
                ss.str("");
                ss.clear();
                ss << *ab;
                auto atms = SymEngine::symbol(ss.str());
                try {
                    auto solns = SymEngine::solve(constraint, atms)->get_args();
                    if (solns.size() == 1) {
                        int idx = -1;
                        for (int ix_ = 0; ix_ < numCoords; ix_++) {
                            if (SymEngine::eq(*coordinates.get(ix_, 0), *atms)) {
                                idx = ix_;
                                break;
                            }
                        }
                        SymEngine::DenseMatrix lo(1, linear_op.ncols());
                        zeros(lo);
                        lo.set(0, idx, SymEngine::integer(1));
                        linear_op.col_join(lo);
                        SymEngine::DenseMatrix no(1, 1);
                        auto negSol = SymEngine::neg((solns[0]));
                        no.set(0, 0, negSol);
                        nonlinear_op.col_join(no);
                        constraint = SymEngine::add(atms, negSol);
                    }
                } catch (std::exception &ee) {
                    logWarn(*constraint, "\t\t", *atms);
                    initial_constraints.push_back(constraint);
                }
            } else {
                initial_constraints.push_back(constraint);
            }
            //Reduce the order of the constraints where possible
            SymEngine::vec_basic partials;
            for (int ix_ = 0; ix_ < numCoords; ix_++) {
                auto c = coordinates.get(ix_, 0);
                ss.str("");
                ss.clear();
                ss << *c;
                auto cs = SymEngine::symbol(ss.str());
                auto cdiff = constraint->diff(cs);
                partials.push_back(cdiff);
            }
            bool reducible = true;
            for (int ix_ = 0; ix_ < offset; ix_++) {
                if (!SymEngine::eq(*partials[ix_], *SymEngine::zero)) {
                    reducible = false;
                    break;
                }
            }
            if (reducible) {
                //If reducible, create rows in linearop+coordinates or nonlinear operators
                //to include the reduced order constraints
                SymEngine::RCP<const SymEngine::Basic> factor = SymEngine::zero;
                SymEngine::RCP<const SymEngine::Basic> nlin = SymEngine::zero;
                std::map<int, SymEngine::RCP<const SymEngine::Basic>> lin_dict;
                for (int ix_ = offset; ix_ < (offset + ss_size); ix_++) {
                    auto coeff = partials[ix_];
                    if (SymEngine::eq(*factor, *SymEngine::zero) && !SymEngine::eq(*coeff, *SymEngine::zero)) {
                        factor = SymEngine::div(SymEngine::integer(1), coeff);
                        lin_dict[ix_] = SymEngine::integer(1);
                    } else if (SymEngine::eq(*factor, *SymEngine::zero) && !SymEngine::eq(*coeff, *SymEngine::zero)) {
                        auto divV = SymEngine::div(coeff, factor);
                        auto new_coeff = SymEngine::evalf(*divV, 53);
                        if (SymEngine::is_a<SymEngine::Integer>(*new_coeff) || SymEngine::is_a<SymEngine::Rational>(*new_coeff)) {
                            lin_dict[ix_] = new_coeff;
                        } else {
                            nlin = SymEngine::add(nlin, SymEngine::mul(new_coeff, coordinates.get(ix_, 0)));
                        }
                    }
                }
                for (int ix_ = (offset + ss_size); ix_ < numCoords; ix_++) {
                    auto coeff = partials[ix_];
                    if (!SymEngine::eq(*coeff, *SymEngine::zero)) {
                        auto cv = coordinates.get(offset + ss_size + ix_, 0);
                        ss.str("");
                        ss.clear();
                        ss << "dot_" << *cv;
                        auto dvc = SymEngine::parse(ss.str());
                        auto coordVec = coordinates.as_vec_basic();
                        auto dc_idx = numCoords;
                        if (std::find(coordVec.begin(), coordVec.end() - cs_size, dvc) == coordVec.end()) {
                            SymEngine::DenseMatrix lo(1, 1);
                            lo.set(0, 0, dvc);
                            coordinates.row_join(lo);
                            cs_size++;
                            numCoords++;
                            n++;
                            SymEngine::DenseMatrix o(linear_op.nrows(), 1);
                            linear_op.row_join(o);
                        } else {
                            dc_idx = std::find(coordVec.begin(), coordVec.end() - cs_size, dvc) - coordVec.begin();
                        }
                        auto eqn = SymEngine::div(coeff, factor);
                        if (SymEngine::is_a<SymEngine::Integer>(*eqn) || SymEngine::is_a<SymEngine::Rational>(*eqn)) {
                            lin_dict[dc_idx] = eqn;
                        } else {
                            nlin = SymEngine::add(nlin, SymEngine::mul(eqn, dvc));
                        }
                    }
                }

                SymEngine::DenseMatrix lo(1, linear_op.ncols());
                zeros(lo);
                for (auto &lin : lin_dict) {
                    lo.set(0, lin.first, lin.second);
                }
                linear_op.col_join(lo);
                SymEngine::DenseMatrix no(1, 1);
                no.set(0, 0, nlin);
                nonlinear_op.col_join(no);

            } else {
                initial_constraints.push_back(constraint);
            }
        }
    }

    //Finally convert the current form of linear and nonlinear operators into smith normal form
    auto snf = getSmithNormalForm(linear_op, nonlinear_op);

    if (initial_constraints.size() > 0)
        snf.constraints.insert(snf.constraints.end(), initial_constraints.begin(), initial_constraints.end());

    return snf;
}

std::tuple<std::string,std::string,char> getDimensions(SymEngine::RCP<const SymEngine::Basic> &eq,
                           std::unordered_map<std::string, std::tuple<std::string, std::string, char>> &dimensions, char vtype = 'u')
{
    if (!SymEngine::is_a<SymEngine::Add>(*SymEngine::expand(eq))) {
        //Check for the condition where the denominator is a sum of terms
        auto ieq = SymEngine::expand(SymEngine::div(SymEngine::one,eq));
        
        if (SymEngine::is_a<SymEngine::Add>(*SymEngine::expand(ieq))) {
            //Find the largest term
            auto terms = SymEngine::expand(ieq)->get_args();
            SymEngine::RCP<const SymEngine::Basic> maxterm;
            size_t tc = 0;
            for (auto t : terms) {
                auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*t);
                if(atoms.size()>tc){
                    tc = atoms.size();
                    maxterm = t;
                }
            }
            //Handle it
            auto neq = SymEngine::div(SymEngine::one,maxterm);
            return getDimensions(neq, dimensions);
        }

        std::ostringstream ss;
        auto equation = SymEngine::eliminateExpLog(eq);
        auto atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*equation);

        SymEngine::map_basic_basic subs;
        for (auto c : atoms) {
            ss.str("");
            ss.clear();
            ss << *c;

            auto dim = std::get<0>(dimensions[ss.str()]);
            if (dim != "") { //Handle dimensionless entities
                try {
                    subs[c] = SymEngine::parse(dim);
                } catch (std::exception &e) {
                    logWarn("Failed to find dimensions for ", ss.str(), dim);
                }
            } else {
                subs[c] = SymEngine::parse("1.0");
            }
        }

        auto res = SymEngine::expand(equation->subs(subs));
        ss.str("");
        ss.clear();
        ss << *res;
        std::string expString = ss.str();

        //Units doesnt seem to handle expressions that start with a negative sign
        //Hack to handle this
        if (expString.c_str()[0] == '-') {
            expString = expString.substr(1);
        }
        auto measure = units::measurement_from_string(expString);
        auto un = measure.units();
        auto mult = un.multiplier();
        auto baseU = un.base_units();
        auto ustring = units::to_string(units::precise_unit(baseU, mult));
        return std::make_tuple(ustring, "0", vtype);
    } else {
        auto terms = SymEngine::expand(eq)->get_args();
        for (auto t : terms) {
            auto c = getDimensions(t, dimensions);
            if (std::get<0>(c) != "") //Return on the first result, expects all terms to be dimensionally correct
                return c;
            else
                logCritical("Dimension could not be deduced for ", *t);
        }
        return std::make_tuple("", "0", vtype);
    }
}

} // namespace BG

