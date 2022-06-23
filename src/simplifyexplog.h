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

#include <symengine/assumptions.h>
#include <symengine/basic.h>
#include <symengine/refine.h>
#include <symengine/visitor.h>

namespace SymEngine {

RCP<const Basic> simplifyExpLog(const RCP<const Basic> &x,
                             const Assumptions *assumptions = nullptr);

RCP<const Basic> eliminateExpLog(const RCP<const Basic> &x,
                             const Assumptions *assumptions = nullptr);                             

class SimplifyExpLogVisitor: public BaseVisitor<SimplifyExpLogVisitor, TransformVisitor>
{
private:
    const Assumptions *assumptions_;

    std::pair<RCP<const Basic>, RCP<const Basic>>
    simplify_pow(const RCP<const Basic> &e, const RCP<const Basic> &b)
    {
        if (is_a<Csc>(*b) and eq(*e, *minus_one)) {
            // csc(expr) ** -1 = sin(expr)
            return std::make_pair(
                one, sin(down_cast<const OneArgFunction &>(*b).get_arg()));
        } else if (is_a<Sec>(*b) and eq(*e, *minus_one)) {
            // sec(expr) ** -1 = cos(expr)
            return std::make_pair(
                one, cos(down_cast<const OneArgFunction &>(*b).get_arg()));
        } else if (is_a<Cot>(*b) and eq(*e, *minus_one)) {
            // cot(expr) ** -1 = tan(expr)
            return std::make_pair(
                one, tan(down_cast<const OneArgFunction &>(*b).get_arg()));
        } else {
            return std::make_pair(e, b);
        }
    }

public:
    using TransformVisitor::bvisit;

    SimplifyExpLogVisitor(const Assumptions *assumptions)
        : BaseVisitor<SimplifyExpLogVisitor, TransformVisitor>()
        , assumptions_(assumptions)
    {
    }

    void bvisit(const Mul &x){
        map_basic_basic map;
        for (const auto &p : x.get_dict()) {
            auto base = apply(p.first);
            //Handle exp
            if(eq(*base,*E)){
                if(!is_a<Log>(*p.second)){
                    //Clean up the scaling for example do the divisions
                    auto ss = SymEngine::expand(p.second);
                    //Apply the exp and reduce log terms, and create Exp for the other terms and create a multiplication
                    //If it is a summation term
                    for(const auto &n: ss->get_args()){
                        bool notReducible = false;
                        if(is_a<Log>(*n)){ //Positive term
                            auto argm = down_cast<const OneArgFunction &>(*n).get_arg(); 
                            Mul::dict_add_term(map, one, argm);                 
                        }else if(is_a<Mul>(*n)){ //Potential negative term
                            auto mterms = n->get_args();
                            vec_basic tmap;
                            SymEngine::RCP<const SymEngine::Basic> lterm;
                            for(const auto &m : mterms){
                                if(is_a<Log>(*m)){
                                    lterm = m;
                                }else{
                                    tmap.push_back(m);
                                }
                            }
                            const auto coeff = mul(tmap);
                            if(!lterm.is_null()){
                                if(eq(*coeff,*minus_one)){
                                    auto argm = down_cast<const OneArgFunction &>(*lterm).get_arg();                  
                                    Mul::dict_add_term(map, minus_one, argm);                 
                                }else{
                                    auto argm = down_cast<const OneArgFunction &>(*lterm).get_arg();                  
                                    Mul::dict_add_term(map, coeff, argm);  
                                }
                            }else{
                                notReducible = true;
                            }
                        }
                        if(notReducible){
                            auto newpair = simplify_pow(n, base);
                            Mul::dict_add_term(map, newpair.first, newpair.second);                         
                        }
                    }
                }else{
                    auto argm = down_cast<const OneArgFunction &>(*p.second).get_arg(); 
                    Mul::dict_add_term(map, one, argm);                      
                }
            }else{
                auto newpair = simplify_pow(p.second, base);
                Mul::dict_add_term(map, newpair.first, newpair.second);  
            }
        }
        result_ = Mul::from_dict(x.get_coef(), std::move(map));
    }

    void bvisit(const Pow &x)
    {
        auto e = apply(x.get_exp());
        auto base = apply(x.get_base());
        auto pair = simplify_pow(e, base);
        result_ = pow(pair.second, pair.first);
    };

    void bvisit(const OneArgFunction &x)
    {
        auto farg = x.get_arg();
        auto newarg = apply(farg);
        result_ = x.create(newarg);
    };
};


class ElminateExpLogVisitor: public BaseVisitor<ElminateExpLogVisitor, TransformVisitor>
{
private:
    const Assumptions *assumptions_;

    std::pair<RCP<const Basic>, RCP<const Basic>>
    simplify_pow(const RCP<const Basic> &e, const RCP<const Basic> &b)
    {
        if (is_a<Csc>(*b) and eq(*e, *minus_one)) {
            // csc(expr) ** -1 = sin(expr)
            return std::make_pair(
                one, sin(down_cast<const OneArgFunction &>(*b).get_arg()));
        } else if (is_a<Sec>(*b) and eq(*e, *minus_one)) {
            // sec(expr) ** -1 = cos(expr)
            return std::make_pair(
                one, cos(down_cast<const OneArgFunction &>(*b).get_arg()));
        } else if (is_a<Cot>(*b) and eq(*e, *minus_one)) {
            // cot(expr) ** -1 = tan(expr)
            return std::make_pair(
                one, tan(down_cast<const OneArgFunction &>(*b).get_arg()));
        } else {
            return std::make_pair(e, b);
        }
    }

public:
    using TransformVisitor::bvisit;

    ElminateExpLogVisitor(const Assumptions *assumptions)
        : BaseVisitor<ElminateExpLogVisitor, TransformVisitor>()
        , assumptions_(assumptions)
    {
    }

    void bvisit(const Mul &x)
    {
        map_basic_basic map;
        for (const auto &p : x.get_dict()) {
            auto base = apply(p.first);
            // //Set exp to be dimensionless
            if (eq(*base, *E)) {
                Mul::dict_add_term(map, one, one);
            } else {
                auto newpair = simplify_pow(p.second, base);
                Mul::dict_add_term(map, newpair.first, newpair.second);
            }
        }
        result_ = Mul::from_dict(x.get_coef(), std::move(map));
    }

    void bvisit(const Pow &x)
    {
        auto e = apply(x.get_exp());
        auto base = apply(x.get_base());
        auto pair = simplify_pow(e, base);
        result_ = pow(pair.second, pair.first);
    };

    void bvisit(const OneArgFunction &x)
    { //Set functions to be dimensionless
        result_ = one;
    };
};

RCP<const Basic> simplifyExpLog(const RCP<const Basic> &x,
                                const Assumptions *assumptions)
{
    auto expr = refine(x, assumptions);
    SimplifyExpLogVisitor b(assumptions);
    return b.apply(expr);
}

RCP<const Basic> eliminateExpLog(const RCP<const Basic> &x,
                             const Assumptions *assumptions){
        auto expr = refine(x, assumptions);
        ElminateExpLogVisitor b(assumptions);
        return simplify(b.apply(expr));                              
}

} // namespace SymEngine
