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
#include <units.hpp>
#include <iostream>
#include <vector>
#include <sstream>

int main(int argc, char *argv[])
{
    auto unit = "10^3 kg*m/s^2";
    auto target = "1 S";
    auto target2 = "N";
    units::precise_unit compare = units::unit_from_string(unit);
    units::precise_unit targetUnit = units::unit_from_string(target);
    std::cout << unit <<" "<< units::to_string(compare)<<std::endl;
    std::cout << target <<" "<< units::to_string(targetUnit)<<" compare "<<(compare.base_units()!=targetUnit.base_units())<<std::endl;
    units::precise_unit targetUnit2 = units::unit_from_string(target2);
    std::cout << target2 <<" "<< units::to_string(targetUnit2)<<" compare "<<(compare!=targetUnit2)<<std::endl;

    std::string av = "microNm";
    units::precise_unit uv = units::unit_from_string(av);
    auto lunit = units::to_string(units::unit_from_string(av));
    std::cout<<av<<"\t"<<lunit<<std::endl;
        auto mult = uv.multiplier();
        auto baseU = uv.base_units();
        std::vector<std::string> unitl;
        std::string multiplier = "";
        if(mult!=1.0){
            multiplier = "multiplier=\""+std::to_string(mult)+"\"";
        }
        if(baseU.meter()!=0){
            std::string mx = multiplier;
            if(baseU.meter()<0){
                mx = "";
            }
            if(baseU.meter()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.meter())+"\" units=\"metre\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"metre\" />");
            }
            if(baseU.meter()>0){
                multiplier = "";
            }            
        }
        if(baseU.kg()!=0){
            std::string mx = multiplier;
            if(baseU.kg()<0){
                mx = "";
            }            
            if(baseU.kg()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.kg())+"\" units=\"kilogram\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"kilogram\" />");
            }
            if(baseU.kg()>0){
                multiplier = "";
            }             
        }
        if(baseU.second()!=0){
            std::string mx = multiplier;
            if(baseU.second()<0){
                mx = "";
            }            
            if(baseU.second()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.second())+"\" units=\"second\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"second\" />");
            }
            if(baseU.second()>0){
                multiplier = "";
            }             
        }
        if(baseU.ampere()!=0){
            std::string mx = multiplier;
            if(baseU.ampere()<0){
                mx = "";
            }            
            if(baseU.ampere()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.ampere())+"\" units=\"ampere\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"ampere\" />");
            }
            if(baseU.ampere()>0){
                multiplier = "";
            }             
        }
        if(baseU.kelvin()!=0){
            std::string mx = multiplier;
            if(baseU.kelvin()<0){
                mx = "";
            }            
            if(baseU.kelvin()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.kelvin())+"\" units=\"kelvin\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"kelvin\" />");
            }
            if(baseU.kelvin()>0){
                multiplier = "";
            }             
        }
        if(baseU.mole()!=0){
            std::string mx = multiplier;
            if(baseU.mole()<0){
                mx = "";
            }            
            if(baseU.mole()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.mole())+"\" units=\"mole\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"mole\" />");
            }
            if(baseU.mole()>0){
                multiplier = "";
            }             
        }        
        if(baseU.candela()!=0){
            std::string mx = multiplier;
            if(baseU.candela()<0){
                mx = "";
            }            
            if(baseU.candela()!=1){
                unitl.push_back("<unit "+mx+" exponent=\""+std::to_string(baseU.candela())+"\" units=\"candela\" />");
            }else{
                unitl.push_back("<unit "+mx+" units=\"candela\" />");
            }
            if(baseU.candela()>0){
                multiplier = "";
            }             
        } 

        std::ostringstream ss;
        ss<<"<units name=\""<<av<<"\">"<<std::endl;
        for(auto c: unitl){
            ss<<"\t"<<c<<std::endl;
        }  
        ss<<"</units>";
        std::string res = ss.str();
        std::cout<<res<<std::endl;    
}