
/**
* @brief Create a Capacitance instance, Linear Capacitor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<Capacitance> createCapacitor(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        
/**
* @brief Create a Inductance instance, Linear Inductor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<Inductance> createInductor(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        
/**
* @brief Create a Resistance instance, Linear Resistor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<Resistance> createResistor(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        
/**
* @brief Create a PotentialSource instance, Constant Voltage Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<PotentialSource> createConstantVoltageSource(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        
/**
* @brief Create a FlowSource instance, Constant Current Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<FlowSource> createConstantCurrentSource(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        
