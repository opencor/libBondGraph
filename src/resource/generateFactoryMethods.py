'''
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
'''
from __future__ import print_function
import json, os, platform
from datetime import datetime
capacitanceFriends = []
resistanceFriends = []
inductanceFriends = []
effortSourceFriends = []
flowSourceFriends = []
bondDimensions = []
createMethods = dict()
domainElementMetadata = dict()
#For standard elements metadata is predefined
domainElementMetadata['Biochemical'] = dict()
domainElementMetadata['Biochemical']['Reaction'] = {
                                    "name":"Reaction",
                                    "description": "Biochemical Reaction",
                                    "class": "passive",
                                    "shortname":"Rx",
                                    "ports":{
                                        "in":{"dir":"in"},
                                        "out":{"dir":"out"}
                                    },                                    
                                    "parameters":{
                                        "r":{"description":"Reaction Rate","dimension":"mol/s","value":0.0},
                                        "R":{"description":"Universal Gas Constant","dimension":"J/K/mol","value":8.31446261815324},
                                        "T":{"description": "Temperature","dimension":"K","value":310.14}
                                    },
                                    "constitutive_relations":[
                                        "f_0 + f_1",
                                        "f_0 - r*(exp(e_0/(R*T)) - exp(e_1/(R*T))"
                                    ]            
                                }
domainElementMetadata['Biochemical']['Concentration'] = {
                                    "name":"Concentration",
                                    "description":"Concentration of Chemical Species",
                                    "class": "passive",
                                    "shortname":"N",
                                    "ports":{
                                        "0":"1"
                                    },                                    
                                    "parameters":{
                                        "k":{"description": "Biochemical Constant; exp(mu_0/RT)/V","dimension":"1/mol","value":1.0},
                                        "R":{"description":"Universal Gas Constant","dimension":"J/K/mol","value":8.31446261815324},
                                        "T":{"description": "Temperature","dimension":"K","value":310.14}
                                    },
                                    "states":{
                                        "a_0":{"description":"Molar Quantity","dimension":"mol","value":1.0}
                                    },
                                    "constitutive_relations":[
                                        "e_0 - R*T*log(k*a_0)",
                                        "f_0 - da_0"
                                    ]            
                                }

domainElementMetadata['Biochemical']['Stoichiometry'] = {
                                    "name":"Stoichiometry",
                                    "description":"Stoichiometry",
                                    "class": "junction",
                                    "shortname":"Y",
                                    "ports":{
                                        "0": {
                                            "dir": "in"
                                        },
                                        "1": {
                                            "dir": "out"
                                        }
                                    },                                    
                                    "parameters":{
                                        "r0":{"description": "Affinity","dimension":"","value":-1.0},
                                        "r1":{"description":"Chemical Power","dimension":"","value":1.0},
                                    },
                                    "constitutive_relations":[
                                        "r0*e_0 + r1*e_1",
                                        "f_0/r0 - f_1/r1"
                                    ]            
                                }
domainElementMetadata['Biochemical']['Chemostat'] = {
          "name": "Chemostat",
          "description": "Constant chemical potential",
          "shortname": "SS",
          "variableprefix": "Ss",
          "class": "passive",
            "ports":{
                "0":"1"
            },           
          "parameters": {
            "mu": {
              "description": "Potential",
              "dimension": "J/mol",
              "value": 1.0
            }
          }
        }

domainElementMetadata['Biochemical']['Flowstat'] = {
          "name": "Flowstat",
          "description": "Constant chemical flux",
          "shortname": "Sf",
          "variableprefix": "Sf",
          "class": "passive",
            "ports":{
                "0":"1"
            },           
          "parameters": {
            "f": {
              "description": "Chemical Flux",
              "dimension": "mol/s",
              "value": 1.0
            }
          }
        }
                            

domainElementMetadata['Generic'] = dict()
domainElementMetadata['Generic']['Transformer'] = {
                                        "name": "Transformer",
                                        "description": "Linear Transformer",
                                        "shortname": "TF",
                                        "variableprefix": "TF",
                                        "class": "junction",
                                        "ports": {
                                                    "in": {
                                                    "effort": {
                                                        "name": "e_",
                                                        "unit": "Volt"
                                                    },
                                                    "flow": {
                                                        "name": "f_",
                                                        "unit": "Ampere"
                                                    },
                                                    "dir": "in"
                                                    },
                                                    "out": {
                                                    "effort": {
                                                        "name": "e_",
                                                        "unit": "Volt"
                                                    },
                                                    "flow": {
                                                        "name": "f_",
                                                        "unit": "Ampere"
                                                    },
                                                    "dir": "out"
                                                    }
                                                },                                        
                                        "parameters": {
                                            "n": {
                                            "description": "Ratio",
                                            "value": 1
                                            }
                                        },
                                        "constitutive_relations": [
                                            "e_1 - n * e_0",
                                            "f_0 + n * f_1"
                                        ]                                        
                                    }
domainElementMetadata['Generic']['Gyrator'] =  {
                                                "name": "Gyrator",
                                                "description": "Linear Gyrator",
                                                "shortname": "GY",
                                                "variableprefix": "GY",
                                                "class": "junction",
                                                "ports": {
                                                            "in": {
                                                            "effort": {
                                                                "name": "e_",
                                                                "unit": "Volt"
                                                            },
                                                            "flow": {
                                                                "name": "f_",
                                                                "unit": "Ampere"
                                                            },
                                                            "dir": "in"
                                                            },
                                                            "out": {
                                                            "effort": {
                                                                "name": "e_",
                                                                "unit": "Volt"
                                                            },
                                                            "flow": {
                                                                "name": "f_",
                                                                "unit": "Ampere"
                                                            },
                                                            "dir": "out"
                                                            }
                                                        },                                                
                                                "parameters": {
                                                    "n": {
                                                    "description": "Ratio",
                                                    "value": 1
                                                    }
                                                },
                                                "constitutive_relations": [
                                                    "e_1 + n*f_0",
                                                    "e_0 - n*f_1"
                                                ]
                                            }      
domainElementMetadata['Generic']['OneJunction'] =  {
                                                "name": "OneJunction",
                                                "description": "Equal Flow Junction",
                                                "shortname": "1",
                                                "variableprefix": "One",
                                                "class": "junction",
                                                "ports": {
                                                    "limit": 100
                                                },
                                                "equate": [
                                                    "f_"
                                                ]
                                            }       
domainElementMetadata['Generic']['ZeroJunction'] =  {
                                                "name": "ZeroJunction",
                                                "description": "Equal Effort Junction",
                                                "shortname": "0",
                                                "variableprefix": "Zero",
                                                "class": "junction",
                                                "ports": {
                                                    "limit": 100
                                                },
                                                "equate": [
                                                    "e_"
                                                ]
                                            }                                                                                                             

friendsfile = f'''#ifndef _FRIENDS_H__
#define _FRIENDS_H__
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

#include "export.h"
#include "thirdparty/json.hpp"

/*
  Generated by 'resource/generateFactorMethods.py' 
  OS : {platform.platform()}
  Release : {platform.release()}
  Version : {platform.version()}
  InstructionSet : {platform.machine()}
*/

namespace BG {{
/**
* Class friend function declarations, macros are called in ElementBase.h, Elements.h and bondgraph.h
*/

/**
Source JSON used to define domains and factoryMethods
//Generated on {datetime.now().strftime("%d %B, %Y %H:%M:%S")}

'''

hfile = '''/*******************************************************************************

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
'''
cppfile = '''/*******************************************************************************

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
'''

def snakeCase(txt_):
    txt = txt_.strip().replace(" ","_")
    if "_" in txt:
        return  ''.join(word.title() for word in txt.split('_'))
    return txt

def generateCapacitance(desc,domain):
    global hfile, cppfile

    #Sanity checks
    try:
        states = desc['states']
        if len(states) != 1:
            print("A capacitance desciption has incorrect number of states. Only single state variable relations are supported.")
            return

        relations = desc['constitutive_relations']
        if len(relations) !=2:
            print(f"A capacitance desciption has incorrect number of constitutive relations. Two equations are required, {len(relations)} provided.")
            return
        stateName = list(states.keys())[0]
        states = states[stateName]
        if not stateName.endswith('_'):
            if '_' in stateName:
                stateName = stateName[0:stateName.index('_')+1]
            else:
                stateName = stateName+"_"

        parameters = desc['parameters']
        name = snakeCase(desc['name'])
        description = desc['description']
        
        decl = f'''
/**
* @brief Create a Capacitance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        '''
        capacitanceFriends.append(f'friend RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a Capacitance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy){{
    auto ptr_ = __createCapacitance_do_not_call_outside_of_lib(proxy);
    RCPLIB::RCP<Capacitance> ptr = RCPLIB::rcp_dynamic_cast<Capacitance>(ptr_);
    ptr->numStates = 1;
    ptr->setDomain("{domain}");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("{stateName}", "{states['value']}"));
    ptr->mParameter.push_back(state);
    state->units = units::to_string(units::unit_from_string("{states['dimension']}"));'''
        for i,param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
            decl += f'''
    RCPLIB::RCP<Value> parameter{i} = RCPLIB::rcp(new Value("{param}", "{pn['value']}"));
    ptr->mParameter.push_back(parameter{i});
    parameter{i}->units = units::to_string(units::unit_from_string("{pn['dimension']}"));'''
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
                    decl += f'''    
    parameter{i}->universalConstant = {isUC};'''
        for i,ceq in enumerate(relations):
            decl += f'''
    ptr->constitutiveEq.push_back("{ceq}");
    ptr->constitutiveEqIndex.push_back(-1);  '''
        decl += f'''
    ptr->setName("{name}");
    return ptr;    
}}  

        '''
        cppfile += decl
        createMethods[domain]["Capacitance"] = f'create{name}'
    except:
        print(f"Error occured while processing capacitance codegen for {domain} domain ")

def generateInductance(desc,domain):
    global hfile, cppfile

    #Sanity checks
    try:
        states = desc['states']
        if len(states) != 1:
            print("A inductance description has incorrect number of states. Only single state variable relations are supported.")
            return

        relations = desc['constitutive_relations']
        if len(relations) !=2:
            print(f"A inductance description has incorrect number of constitutive relations. Two equations are required, {len(relations)} provided.")
            return
        stateName = list(states.keys())[0]
        states = states[stateName]
        if not stateName.endswith('_'):
            if '_' in stateName:
                stateName = stateName[0:stateName.index('_')+1]
            else:
                stateName = stateName+"_"

        parameters = desc['parameters']
        name = snakeCase(desc['name'])
        description = desc['description']

        decl = f'''
/**
* @brief Create a Inductance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        '''
        inductanceFriends.append(f'friend RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a Inductance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy){{
    auto ptr_ = __createInductance_do_not_call_outside_of_lib(proxy);
    RCPLIB::RCP<Inductance> ptr = RCPLIB::rcp_dynamic_cast<Inductance>(ptr_);
    ptr->numStates = 1;
    ptr->setDomain("{domain}");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("{stateName}", "{states['value']}"));
    ptr->mParameter.push_back(state);
    state->units = units::to_string(units::unit_from_string("{states['dimension']}"));'''
        for i,param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    RCPLIB::RCP<Value> parameter{i} = RCPLIB::rcp(new Value("{param}", "{pn['value']}"));
    ptr->mParameter.push_back(parameter{i});
    parameter{i}->units = units::to_string(units::unit_from_string("{pn['dimension']}"));'''
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
                    decl += f'''
    parameter{i}->universalConstant = {isUC};'''
        for i,ceq in enumerate(relations):
            decl += f'''
    ptr->constitutiveEq.push_back("{ceq}");
    ptr->constitutiveEqIndex.push_back(-1);  '''
        decl += f'''
    ptr->setName("{name}");
    return ptr;
}}  
        '''
        cppfile += decl
        createMethods[domain]["Inductance"] = f'create{name}'
    except:
        print(f"Error occured while processing inductance codegen for {domain} domain ")

def generateResistance(desc,domain):
    global hfile, cppfile

    #Sanity checks
    try:
        relations = desc['constitutive_relations']
        if len(relations) !=1:
            print(f"A resistance desciption has incorrect number of constitutive relations. One equation is required, {len(relations)} provided.")
            return

        parameters = desc['parameters']
        name = snakeCase(desc['name'])
        description = desc['description']

        decl = f'''
/**
* @brief Create a Resistance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        '''
        resistanceFriends.append(f'friend RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a Resistance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy){{
    auto ptr_ = __createResistance_do_not_call_outside_of_lib(proxy);
    RCPLIB::RCP<Resistance> ptr = RCPLIB::rcp_dynamic_cast<Resistance>(ptr_);
    ptr->numStates = 0;
    ptr->setDomain("{domain}");'''
        for i, param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    RCPLIB::RCP<Value> parameter{i} = RCPLIB::rcp(new Value("{param}", "{pn['value']}"));
    ptr->mParameter.push_back(parameter{i});
    parameter{i}->units = units::to_string(units::unit_from_string("{pn['dimension']}"));'''
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
                    decl += f'''
    parameter{i}->universalConstant = {isUC};'''
        for i,ceq in enumerate(relations):
            decl += f'''
    ptr->constitutiveEq.push_back("{ceq}");
    ptr->constitutiveEqIndex.push_back(-1);  '''
        decl += f'''
    ptr->setName("{name}");
    return ptr;
}}  

        '''
        cppfile += decl
        createMethods[domain]["Resistance"] = f'create{name}'
    except:
        print(f"Error occured while processing Resistance codegen for {domain} domain ")

def generatePotentialSource(desc,domain):
    global hfile, cppfile

    #Sanity checks
    try:
        parameters = desc['parameters']
        name = snakeCase(desc['name'])
        description = desc['description']

        decl = f'''
/**
* @brief Create a PotentialSource instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        '''
        effortSourceFriends.append(f'friend RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a PotentialSource instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy){{
    auto ptr_ = createPotentialSource(proxy);
    RCPLIB::RCP<PotentialSource> ptr = RCPLIB::rcp_dynamic_cast<PotentialSource>(ptr_);
    ptr->setDomain("{domain}");'''
        for i, param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    ptr->mParameter[{i}]->name = "{param}";
    ptr->mParameter[{i}]->value = "{pn['value']}";
    ptr->mParameter[{i}]->units = units::to_string(units::unit_from_string("{pn['dimension']}"));'''
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
                    decl += f'''
    ptr->parameter{i}->universalConstant = {isUC};'''
        decl += f'''
    ptr->setName("{name}");
    return ptr;    
}}  

        '''
        cppfile += decl
        createMethods[domain]["PotentialSource"] = f'create{name}'
    except:
        print(f"Error occured while processing PotentialSource codegen for {domain} domain ")

def generateFlowSource(desc,domain):
    global hfile, cppfile

    #Sanity checks
    try:
        parameters = desc['parameters']
        name = snakeCase(desc['name'])
        description = desc['description']

        decl = f'''
/**
* @brief Create a FlowSource instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
EXPORTED RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy=RCPLIB::null);

        '''
        flowSourceFriends.append(f'friend RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a FlowSource instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<BGElement> create{name}(const RCPLIB::RCP< BGElementData > &proxy){{
    auto ptr_ = createFlowSource(proxy);
    RCPLIB::RCP<FlowSource> ptr = RCPLIB::rcp_dynamic_cast<FlowSource>(ptr_);
    ptr->setDomain("{domain}");'''
        for i, param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    ptr->mParameter[{i}]->name = "{param}";
    ptr->mParameter[{i}]->value = "{pn['value']}";
    ptr->mParameter[{i}]->units = units::to_string(units::unit_from_string("{pn['dimension']}"));'''
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
                    decl += f'''
    ptr->parameter{i}->universalConstant = {isUC};'''

        decl += f'''
    ptr->setName("{name}");
    return ptr;    
}}  

        '''
        cppfile += decl
        createMethods[domain]["FlowSource"] = f'create{name}'
    except:
        print(f"Error occured while processing FlowSource codegen for {domain} domain ")


def generateStatementsForDomain(relations,domain):
    createMethods[domain] = dict()
    if "Capacitance" in relations:
        generateCapacitance(relations["Capacitance"],domain)
    if "Inductance" in relations:
        generateInductance(relations["Inductance"], domain)
    if "Resistance" in relations:
        generateResistance(relations["Resistance"], domain)
    if "PotentialSource" in relations:
        generatePotentialSource(relations["PotentialSource"], domain)
    if "FlowSource" in relations:
        generateFlowSource(relations["FlowSource"], domain)
    if "BondDimensions" in relations:
        br = [domain]
        desc = relations['BondDimensions']
        br.append(desc['effort'])
        br.append(desc['flow'])
        bondDimensions.append(br)

def generateMetaData(relations,domain):
    if not domain in domainElementMetadata:
        domainElementMetadata[domain] = dict()
    if "Capacitance" in relations:
        domainElementMetadata[domain]["Capacitance"] = relations["Capacitance"]
    if "Inductance" in relations:
        domainElementMetadata[domain]["Inductance"] = relations["Inductance"]
    if "Resistance" in relations:
        domainElementMetadata[domain]["Resistance"] = relations["Resistance"]
    if "PotentialSource" in relations:
        domainElementMetadata[domain]["PotentialSource"] = relations["PotentialSource"]
    if "FlowSource" in relations:
        domainElementMetadata[domain]["FlowSource"] = relations["FlowSource"]
    

if __name__ == '__main__':
    spath = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(spath,'constitutiveRelations.json'),'r') as jfile:
        constitutiverelationsAll = json.load(jfile)
        if 'ElementDefinitions' in constitutiverelationsAll:
            constitutiverelations = constitutiverelationsAll['ElementDefinitions']
            for domain,desc in constitutiverelations.items():
                if type(desc) is dict:
                    generateStatementsForDomain(desc, domain)
                    generateMetaData(desc,domain)
            with open(os.path.join(spath,'../factorymethods.h'), 'w') as hf:
                print(hfile,file=hf)
            with open(os.path.join(spath,'../factorymethods.cpp'), 'w') as hc:
                print(cppfile,file=hc)

            cstmts = ['''/*******************************************************************************

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
                      
#include "thirdparty/json.hpp"
//metadatastring is generated by cmake, kept out of namespace as vtk_encode_string does not have a option
//to specify namespace
extern const char *metadatastring;
namespace BG {
                    '''] # for friends.cpp file
            with open(os.path.join(spath,'../friends.h'),'w') as frd:            
                stmts = [friendsfile+json.dumps(constitutiverelations,indent=4, sort_keys=True)]
               
                stmts.append('\n*/')
                if len(capacitanceFriends)>0:
                    res = "#define DEFINE_FRIENDS_OF_CAPACITANCE \\"
                    if len(capacitanceFriends) > 1:
                        for cp in capacitanceFriends[:-1]:
                            res += '\n\t'+cp+'\t\\'
                    res += '\n\t'+capacitanceFriends[-1]
                    stmts.append(res)
                else:
                    stmts.append("#define DEFINE_FRIENDS_OF_CAPACITANCE ")

                if len(inductanceFriends)>0:
                    res = "#define DEFINE_FRIENDS_OF_INDUCTANCE \\"
                    if len(inductanceFriends) > 1:
                        for cp in inductanceFriends[:-1]:
                            res += '\n\t'+cp+'\t\\'
                    res += '\n\t'+inductanceFriends[-1]
                    stmts.append(res)
                else:
                    stmts.append("#define DEFINE_FRIENDS_OF_INDUCTANCE")

                if len(resistanceFriends)>0:
                    res = "#define DEFINE_FRIENDS_OF_RESISTANCE \\"
                    if len(resistanceFriends) > 1:
                        for cp in resistanceFriends[:-1]:
                            res += '\n\t'+cp+'\t\\'
                    res += '\n\t'+resistanceFriends[-1]
                    stmts.append(res)
                else:
                    stmts.append("#define DEFINE_FRIENDS_OF_RESISTANCE")

                if len(effortSourceFriends)>0:
                    res = "#define DEFINE_FRIENDS_OF_EFFORTSOURCE \\"
                    if len(effortSourceFriends) > 1:
                        for cp in effortSourceFriends[:-1]:
                            res += '\n\t'+cp+'\t\\'
                    res += '\n\t'+effortSourceFriends[-1]
                    stmts.append(res)
                else:
                    stmts.append("#define DEFINE_FRIENDS_OF_EFFORTSOURCE")

                if len(flowSourceFriends)>0:
                    res = "#define DEFINE_FRIENDS_OF_FLOWSOURCE \\"
                    if len(flowSourceFriends) > 1:
                        for cp in flowSourceFriends[:-1]:
                            res += '\n\t'+cp+'\t\\'
                    res += '\n\t'+flowSourceFriends[-1]
                    stmts.append(res)
                else:
                    stmts.append("#define DEFINE_FRIENDS_OF_FLOWSOURCE")
                stmts.append("\n\n//Physical domain bond dimensions")
                stmts.append("#define DEFINE_PHYSICAL_DOMAINS \\")
                stmts.append("\tstd::unordered_map<std::string,std::tuple<std::string,std::string> > portDimensions = { \\")
                for bd in bondDimensions[:-1]:
                    stmts.append(f'\t\t{{"{bd[0]}",std::make_tuple("{bd[1]}","{bd[2]}")}}, \\')
                bd = bondDimensions[-1]
                stmts.append(f'\t\t{{"{bd[0]}",std::make_tuple("{bd[1]}","{bd[2]}")}}\\')
                stmts.append("\t};")

                cstmts.append("\n\nnlohmann::json getSupportedPhysicalDomainsAndFactoryMethods() {")
                cstmts.append("\tnlohmann::json result;")
                dlist = '"Biochemical",'

                #Add standard ones for biochemical and generic
                if not 'Biochemical' in createMethods:
                    createMethods['Biochemical'] = dict()
                createMethods['Biochemical']["Reaction"] = "createReaction"
                createMethods['Biochemical']["Concentration"] = "createConcentration"
                createMethods['Biochemical']["Chemostat"] = "createChemostat"
                createMethods['Biochemical']["Flowstat"] = "createFlowstat"

                if not 'Generic' in createMethods:
                    createMethods['Generic'] = dict()
                createMethods['Generic']["Transformer"] = "createTransformer"
                createMethods['Generic']["Gyrator"] = "createGyrator"
                createMethods['Generic']["OneJunction"] = "createOneJunction"
                createMethods['Generic']["ZeroJunction"] = "createZeroJunction"
                #createMethods['Generic']["PotentialSource"] = "createConstantPotentialSource"
                #createMethods['Generic']["FlowSource"] = "createConstantFlowSource"

                if not 'Composition' in createMethods:
                    createMethods['Composition'] = dict()
                createMethods['Composition']["PortHamiltonian"] = "createPortHamiltonian"
                createMethods['Composition']["UserDefined"] = "createUserDefined"

                mlist = []
                for d,mts in createMethods.items():
                    if len(mts) > 0:
                        dlist += f'"{d}",'
                        cstmts.append(f"\t\n\tnlohmann::json mts{d};")
                        for md,mn in mts.items():
                            cstmts.append(f'\tmts{d}["{md}"]="{mn}";')
                            mlist.append(mn)
                        cstmts.append(f'\tresult["{d}"] = mts{d};\n\t')

                cstmts.append(f'\tresult["Domains"] = nlohmann::json::array({{{dlist[:-1]}}});\n')
                cstmts.append('\tresult["ElementDefinitions"]  = nlohmann::json::parse(std::string(metadatastring));');
                cstmts.append("\treturn result;")
                cstmts.append("}\n}")

                stmts.append('\n\n//Reversing strings to optimise string comparisions as most methods start with create')
                stmts.append("#define CALL_FACTORY_METHODS_BY_NAME \\")
                stmts.append("\tstd::string rMethodName(methodName.rbegin(),methodName.rend()); \\")
                stmts.append("\tRCPLIB::RCP<BGElementData> data = RCPLIB::null;  \\")
                stmts.append(f'\tif(rMethodName=="{mlist[0][::-1]}") {{return {mlist[0]}(data);}} \\')
                for m in mlist[1:]:
                    if m!="createPortHamiltonian":
                        stmts.append(f'\telse if(rMethodName=="{m[::-1]}") {{return {m}(data);}} \\')
                    else:
                        #PH requires an argument, this is passed through the function parameter `expr`
                        stmts.append(f'\telse if(rMethodName=="{m[::-1]}") {{nlohmann::json j = nlohmann::json::parse(expr); std::vector<std::string> states = j["states"]; return {m}(j["Hamiltonian"],states,data);}} \\')

                stmts.append("\n}\n#endif") #End header
                print('\n'.join(stmts),file=frd)
            with open(os.path.join(spath,'../friends.cpp'),'w') as frdcpp:
                print('\n'.join(cstmts),file=frdcpp)

            with open(os.path.join(spath,'metadata.json'),'w') as meta:
                print(json.dumps(domainElementMetadata,indent=4),file=meta)
        else:
            print("Incorrect format of constitutive relations file!")
            import sys
            sys.exit(-1)