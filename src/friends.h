#ifndef _FRIENDS_H__
#define _FRIENDS_H__
#include "export.h"
#include "thirdparty/json.hpp"
/**
* Class friend function declarations, macros are called in ElementBase.h, Elements.h and bondgraph.h
*/

/**
Source JSON used to define domains and factoryMethods
//Generated on 29 July, 2021 22:16:04

{
    "Electrical": {
        "BondDimensions": {
            "effort": "V",
            "flow": "A"
        },
        "Capacitance": {
            "constitutive_relations": [
                "q_0 - C * e_0",
                "dot_q_0 - f_0"
            ],
            "description": "Linear Capacitor",
            "name": "Capacitor",
            "parameters": {
                "C": {
                    "description": "Capacitance",
                    "dimension": "Farad",
                    "value": 1.0
                }
            },
            "states": {
                "q_0": {
                    "description": "Electric charge",
                    "dimension": "coulomb",
                    "value": 0.0
                }
            }
        },
        "PotentialSource": {
            "description": "Constant Voltage Source",
            "name": "ConstantVoltageSource",
            "parameters": {
                "u": {
                    "description": "Voltage",
                    "dimension": "V",
                    "value": 1.0
                }
            }
        },
        "FlowSource": {
            "description": "Constant Current Source",
            "name": "ConstantCurrentSource",
            "parameters": {
                "a": {
                    "description": "Current",
                    "dimension": "A",
                    "value": 1.0
                }
            }
        },
        "Inductance": {
            "constitutive_relations": [
                "p_0 - L*f_0",
                "dot_p_0 - e_0"
            ],
            "description": "Linear Inductor",
            "name": "Inductor",
            "parameters": {
                "L": {
                    "description": "Inductance",
                    "dimension": "henry",
                    "value": 1.0
                }
            },
            "states": {
                "p_0": {
                    "description": "Magetic flux",
                    "dimension": "weber",
                    "value": 0.0
                }
            }
        },
        "Resistance": {
            "constitutive_relations": [
                "e_0 - f_0*r"
            ],
            "description": "Linear Resistor",
            "name": "Resistor",
            "parameters": {
                "r": {
                    "description": "Resistance",
                    "dimension": "Ohm",
                    "universalConstant": false,
                    "value": 1.0
                }
            }
        }
    }
}

*/
#define DEFINE_FRIENDS_OF_CAPACITANCE \
	friend RCPLIB::RCP<Capacitance> createCapacitor(const RCPLIB::RCP< ElementImpl > &proxy);
#define DEFINE_FRIENDS_OF_INDUCTANCE \
	friend RCPLIB::RCP<Inductance> createInductor(const RCPLIB::RCP< ElementImpl > &proxy);
#define DEFINE_FRIENDS_OF_RESISTANCE \
	friend RCPLIB::RCP<Resistance> createResistor(const RCPLIB::RCP< ElementImpl > &proxy);
#define DEFINE_FRIENDS_OF_PotentialSource \
	friend RCPLIB::RCP<PotentialSource> createConstantVoltageSource(const RCPLIB::RCP< ElementImpl > &proxy);
#define DEFINE_FRIENDS_OF_FLOWSOURCE \
	friend RCPLIB::RCP<FlowSource> createConstantCurrentSource(const RCPLIB::RCP< ElementImpl > &proxy);


//Physical domain bond dimensions
#define DEFINE_PHYSICAL_DOMAINS \
	std::unordered_map<std::string,std::tuple<std::string,std::string> > portDimensions = { \
		{"Electrical",std::make_tuple("V","A")}\
	};


EXPORTED inline nlohmann::json getSupportedPhysicalDomainsAndFactoryMethods() {
	nlohmann::json result;
	
	nlohmann::json mtsElectrical;
	mtsElectrical["Capacitance"]="createCapacitor";
	mtsElectrical["Inductance"]="createInductor";
	mtsElectrical["Resistance"]="createResistor";
	mtsElectrical["PotentialSource"]="createConstantVoltageSource";
	mtsElectrical["FlowSource"]="createConstantCurrentSource";
	result["Electrical"] = mtsElectrical;
	
	result["Domains"] = nlohmann::json::array({"Biochemical","Electrical"});
	return result;
}


//Reversing strings to optimise string comparisions as most methods start with create
#define CALL_FACTORY_METHODS_BY_NAME \
	std::string rMethodName(methodName.rbegin(),methodName.rend()); \
	RCPLIB::RCP<ElementImpl> data = RCPLIB::null;  \
	if(rMethodName=="roticapaCetaerc") {return createCapacitor(data);} \
	else if(rMethodName=="rotcudnIetaerc") {return createInductor(data);} \
	else if(rMethodName=="rotsiseRetaerc") {return createResistor(data);} \
	else if(rMethodName=="ecruoSegatloVtnatsnoCetaerc") {return createConstantVoltageSource(data);} \
	else if(rMethodName=="ecruoStnerruCtnatsnoCetaerc") {return createConstantCurrentSource(data);} \

#endif
