
/**
* @brief Create a Capacitance instance, Linear Capacitor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<Capacitance> createCapacitor(const RCPLIB::RCP< ElementImpl > &proxy){
    auto ptr = __createCapacitance_do_not_call_outside_of_lib(proxy);
    ptr->numStates = 1;
    ptr->setDomain("Electrical");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("q_", "0.0"));
    ptr->mParameter.push_back(state);
    state->units = units::unit_from_string("coulomb");
    RCPLIB::RCP<Value> parameter0 = RCPLIB::rcp(new Value("C", "1.0"));
    ptr->mParameter.push_back(parameter0);
    parameter0->units = units::unit_from_string("Farad");
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
RCPLIB::RCP<Inductance> createInductor(const RCPLIB::RCP< ElementImpl > &proxy){
    auto ptr = __createInductance_do_not_call_outside_of_lib(proxy);
    ptr->numStates = 1;
    ptr->setDomain("Electrical");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("p_", "0.0"));
    ptr->mParameter.push_back(state);
    state->units = units::unit_from_string("weber");
    RCPLIB::RCP<Value> parameter0 = RCPLIB::rcp(new Value("L", "1.0"));
    ptr->mParameter.push_back(parameter0);
    parameter0->units = units::unit_from_string("henry");
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
RCPLIB::RCP<Resistance> createResistor(const RCPLIB::RCP< ElementImpl > &proxy){
    auto ptr = __createResistance_do_not_call_outside_of_lib(proxy);
    ptr->numStates = 0;
    ptr->setDomain("Electrical");
    RCPLIB::RCP<Value> parameter0 = RCPLIB::rcp(new Value("r", "1.0"));
    ptr->mParameter.push_back(parameter0);
    parameter0->units = units::unit_from_string("Ohm");
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
RCPLIB::RCP<PotentialSource> createConstantVoltageSource(const RCPLIB::RCP< ElementImpl > &proxy){
    auto ptr = createPotentialSource(proxy);
    ptr->setDomain("Electrical");
    ptr->mParameter[0]->name = "u";
    ptr->mParameter[0]->value = "1.0";
    ptr->mParameter[0]->units = units::unit_from_string("V");
    ptr->setName("ConstantVoltageSource");
    return ptr;    
}  

        
/**
* @brief Create a FlowSource instance, Constant Current Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<FlowSource> createConstantCurrentSource(const RCPLIB::RCP< ElementImpl > &proxy){
    auto ptr = createFlowSource(proxy);
    ptr->setDomain("Electrical");
    ptr->mParameter[0]->name = "a";
    ptr->mParameter[0]->value = "1.0";
    ptr->mParameter[0]->units = units::unit_from_string("A");
    ptr->setName("ConstantCurrentSource");
    return ptr;    
}  

        
