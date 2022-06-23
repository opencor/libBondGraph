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
* @brief Create a Capacitance instance, Linear Capacitor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createCapacitor(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = __createCapacitance_do_not_call_outside_of_lib(proxy);
    RCPLIB::RCP<Capacitance> ptr = RCPLIB::rcp_dynamic_cast<Capacitance>(ptr_);
    ptr->numStates = 1;
    ptr->setDomain("Electrical");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("q_", "0"));
    ptr->mParameter.push_back(state);
    state->units = units::to_string(units::unit_from_string("coulomb"));
    RCPLIB::RCP<Value> parameter0 = RCPLIB::rcp(new Value("C", "1"));
    ptr->mParameter.push_back(parameter0);
    parameter0->units = units::to_string(units::unit_from_string("Farad"));
    ptr->constitutiveEq.push_back("q_0 - C * e_0");
    ptr->constitutiveEqIndex.push_back(-1);  
    ptr->constitutiveEq.push_back("dot_q_0 - f_0");
    ptr->constitutiveEqIndex.push_back(-1);  
    ptr->setName("Capacitor");
    return ptr;    
}  

        
/**
* @brief Create a Inductance instance, Linear Inductor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createInductor(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = __createInductance_do_not_call_outside_of_lib(proxy);
    RCPLIB::RCP<Inductance> ptr = RCPLIB::rcp_dynamic_cast<Inductance>(ptr_);
    ptr->numStates = 1;
    ptr->setDomain("Electrical");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("p_", "0"));
    ptr->mParameter.push_back(state);
    state->units = units::to_string(units::unit_from_string("weber"));
    RCPLIB::RCP<Value> parameter0 = RCPLIB::rcp(new Value("L", "1"));
    ptr->mParameter.push_back(parameter0);
    parameter0->units = units::to_string(units::unit_from_string("henry"));
    ptr->constitutiveEq.push_back("p_0 - L*f_0");
    ptr->constitutiveEqIndex.push_back(-1);  
    ptr->constitutiveEq.push_back("dot_p_0 - e_0");
    ptr->constitutiveEqIndex.push_back(-1);  
    ptr->setName("Inductor");
    return ptr;
}  
        
/**
* @brief Create a Resistance instance, Linear Resistor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createResistor(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = __createResistance_do_not_call_outside_of_lib(proxy);
    RCPLIB::RCP<Resistance> ptr = RCPLIB::rcp_dynamic_cast<Resistance>(ptr_);
    ptr->numStates = 0;
    ptr->setDomain("Electrical");
    RCPLIB::RCP<Value> parameter0 = RCPLIB::rcp(new Value("r", "1"));
    ptr->mParameter.push_back(parameter0);
    parameter0->units = units::to_string(units::unit_from_string("Ohm"));
    ptr->constitutiveEq.push_back("e_0 - f_0*r");
    ptr->constitutiveEqIndex.push_back(-1);  
    ptr->setName("Resistor");
    return ptr;
}  

        
/**
* @brief Create a PotentialSource instance, Constant Voltage Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createConstantVoltageSource(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = createPotentialSource(proxy);
    RCPLIB::RCP<PotentialSource> ptr = RCPLIB::rcp_dynamic_cast<PotentialSource>(ptr_);
    ptr->setDomain("Electrical");
    ptr->mParameter[0]->name = "u";
    ptr->mParameter[0]->value = "1";
    ptr->mParameter[0]->units = units::to_string(units::unit_from_string("V"));
    ptr->setName("ConstantVoltageSource");
    return ptr;    
}  

        
/**
* @brief Create a FlowSource instance, Constant Current Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createConstantCurrentSource(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = createFlowSource(proxy);
    RCPLIB::RCP<FlowSource> ptr = RCPLIB::rcp_dynamic_cast<FlowSource>(ptr_);
    ptr->setDomain("Electrical");
    ptr->mParameter[0]->name = "i";
    ptr->mParameter[0]->value = "1";
    ptr->mParameter[0]->units = units::to_string(units::unit_from_string("A"));
    ptr->setName("ConstantCurrentSource");
    return ptr;    
}  

        
/**
* @brief Create a PotentialSource instance, Constant Potential Source
* Physical domain Generic
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createConstantPotentialSource(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = createPotentialSource(proxy);
    RCPLIB::RCP<PotentialSource> ptr = RCPLIB::rcp_dynamic_cast<PotentialSource>(ptr_);
    ptr->setDomain("Generic");
    ptr->mParameter[0]->name = "u";
    ptr->mParameter[0]->value = "1";
    ptr->mParameter[0]->units = units::to_string(units::unit_from_string("V"));
    ptr->setName("ConstantPotentialSource");
    return ptr;    
}  

        
/**
* @brief Create a FlowSource instance, Constant Flow Source
* Physical domain Generic
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> createConstantFlowSource(const RCPLIB::RCP< BGElementData > &proxy){
    auto ptr_ = createFlowSource(proxy);
    RCPLIB::RCP<FlowSource> ptr = RCPLIB::rcp_dynamic_cast<FlowSource>(ptr_);
    ptr->setDomain("Generic");
    ptr->mParameter[0]->name = "i";
    ptr->mParameter[0]->value = "1";
    ptr->mParameter[0]->units = units::to_string(units::unit_from_string("A"));
    ptr->setName("ConstantFlowSource");
    return ptr;    
}  

        
