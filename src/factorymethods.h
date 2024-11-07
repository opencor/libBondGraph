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
EXPORTED RCPLIB::RCP<BGElement> createCapacitor(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Inductance instance, Linear Inductor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createInductor(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Resistance instance, Linear Resistor
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createResistor(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a PotentialSource instance, Constant Voltage Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantVoltageSource(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a FlowSource instance, Constant Current Source
* Physical domain Electrical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantCurrentSource(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Capacitance instance, Linear Elastic Spring
* Physical domain Mechanical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createLinearSpring(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Inductance instance, Mass
* Physical domain Mechanical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createMass(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Resistance instance, Linear damping
* Physical domain Mechanical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createLinearDamper(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a PotentialSource instance, Constant Point Force
* Physical domain Mechanical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantForce(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a FlowSource instance, Constant Linear Velocity
* Physical domain Mechanical
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantVelocity(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Capacitance instance, Vessal wall compliance
* Physical domain Hydraulic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createFluidCompliance(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Inductance instance, Fluid Inertia
* Physical domain Hydraulic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createFluidInertance(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a Resistance instance, Viscous Resistance
* Physical domain Hydraulic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createViscousResistance(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a PotentialSource instance, Constant Pressure Source
* Physical domain Hydraulic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantPressureSource(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a FlowSource instance, Constant Fluid Flow Source
* Physical domain Hydraulic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantFluidFlowSource(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a PotentialSource instance, Constant Potential Source
* Physical domain Generic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantPotentialSource(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
/**
* @brief Create a FlowSource instance, Constant Flow Source
* Physical domain Generic
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> createConstantFlowSource(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        
