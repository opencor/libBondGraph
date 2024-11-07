
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

*******************************************************************************/#include <symengine/visitor.h>
#include <symengine/basic.h>
#include <symengine/assumptions.h>
#include <symengine/refine.h>

namespace SymEngine
{

RCP<const Basic> simplifyExpLog(const RCP<const Basic> &x,
                          const Assumptions *assumptions = nullptr);


RCP<const Basic> eliminateExpLog(const RCP<const Basic> &x,
                             const Assumptions *assumptions = nullptr);                            

class SimplifyExpLogVisitor : public BaseVisitor<SimplifyExpLogVisitor, TransformVisitor>
{
private:
    const Assumptions *assumptions_;

    std::pair<RCP<const Basic>, RCP<const Basic>>
    simplify_pow(const RCP<const Basic> &e, const RCP<const Basic> &b){
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
        : BaseVisitor<SimplifyExpLogVisitor, TransformVisitor>(),
          assumptions_(assumptions)
    {
    }

    void bvisit(const Mul &x){
        map_basic_basic map;
        for (const auto &p : x.get_dict()) {
            auto base = apply(p.first);
            //Handle exp
            std::cout<<*p.first<<"\t"<<*p.second<<std::endl;
            if(eq(*base,*E)){
                if(!is_a<Log>(*p.second)){
                    //Clean up the scaling for example do the divisions
                    auto ss = SymEngine::expand(p.second);
                    //Apply the exp and reduce log terms, and create Exp for the other terms and create a multiplication
                    //If it is a summation term
                    for(const auto &n: ss->get_args()){
                        bool notReducible = false;
                        std::cout<<"Args \t"<<*n<<std::endl;
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

    void bvisit(const OneArgFunction &x){
        auto farg = x.get_arg();
        auto newarg = apply(farg);
        result_ = x.create(newarg);
    };
};

RCP<const Basic> simplifyExpLog(const RCP<const Basic> &x,
                          const Assumptions *assumptions){
        auto expr = refine(x, assumptions);
        SimplifyExpLogVisitor b(assumptions);
        return b.apply(expr);                              
    }


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

RCP<const Basic> eliminateExpLog(const RCP<const Basic> &x,
                             const Assumptions *assumptions){
        auto expr = refine(x, assumptions);
        ElminateExpLogVisitor b(assumptions);
        return b.apply(expr);                              
}


}

using namespace SymEngine;

#include <iostream>
#include <symengine/parser.h>

int main(int argc, char *argv[])
{
    auto x = symbol("x");

    // The following tests visits Mul
    {
    //auto expr = SymEngine::parse("R*T*log(k_10*a_10)");
    //std::cout<<*expr<<"\t"<<*simplifyExp(expr)<<std::endl;
    }
    {
        //auto expr = SymEngine::parse("exp((1/(R*T))*(R*T*log(k_10*a_10) + R*T*log(k_11*a_11)))*r_12");
        //std::cout<<*expr<<"\t"<<*simplifyExpLog(expr)<<std::endl;
        auto expr = SymEngine::parse("-exp(-(-R*T*log(a_5*k_5) - R*T*log(a_6*k_6))/(R*T))*r_19 - exp(-(-R*T*log(a_5*k_5) - R*T*log(k_13*a_13))/(R*T))*r_18 + exp(log(a_0*k_0))*r_19 + exp(log(a_4*k_4))*r_18");
        std::cout<<*expr<<"\t"<<*simplifyExpLog(expr)<<std::endl;

        //auto expr = SymEngine::parse("-R*T*log(a_5*k_5)");
        //std::cout<<*expr<<"\t"<<*eliminateExpLog(expr)<<std::endl;
    }

}