from __future__ import print_function
import json, os
from datetime import datetime
capacitanceFriends = []
resistanceFriends = []
inductanceFriends = []
PotentialSourceFriends = []
flowSourceFriends = []
bondDimensions = []
createMethods = dict()

friendsfile = f'''#ifndef _FRIENDS_H__
#define _FRIENDS_H__
#include "export.h"
#include "thirdparty/json.hpp"
/**
* Class friend function declarations, macros are called in ElementBase.h, Elements.h and bondgraph.h
*/

/**
Source JSON used to define domains and factoryMethods
//Generated on {datetime.now().strftime("%d %B, %Y %H:%M:%S")}

'''

hfile = ''
cppfile = ''

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
EXPORTED RCPLIB::RCP<Capacitance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        '''
        capacitanceFriends.append(f'friend RCPLIB::RCP<Capacitance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a Capacitance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<Capacitance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy){{
    auto ptr = __createCapacitance_do_not_call_outside_of_lib(proxy);
    ptr->numStates = 1;
    ptr->setDomain("{domain}");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("{stateName}", "{states['value']}"));
    ptr->mParameter.push_back(state);
    state->units = units::unit_from_string("{states['dimension']}");'''
        for i,param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            if 'universalConstant' in pn:
                if pn['universalConstant']:
                    isUC = "true"
            decl += f'''
    RCPLIB::RCP<Value> parameter{i} = RCPLIB::rcp(new Value("{param}", "{pn['value']}"));
    ptr->mParameter.push_back(parameter{i});
    parameter{i}->units = units::unit_from_string("{pn['dimension']}");'''
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
EXPORTED RCPLIB::RCP<Inductance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        '''
        inductanceFriends.append(f'friend RCPLIB::RCP<Inductance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a Inductance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<Inductance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy){{
    auto ptr = __createInductance_do_not_call_outside_of_lib(proxy);
    ptr->numStates = 1;
    ptr->setDomain("{domain}");
    RCPLIB::RCP<Value> state = RCPLIB::rcp(new Value("{stateName}", "{states['value']}"));
    ptr->mParameter.push_back(state);
    state->units = units::unit_from_string("{states['dimension']}");'''
        for i,param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    RCPLIB::RCP<Value> parameter{i} = RCPLIB::rcp(new Value("{param}", "{pn['value']}"));
    ptr->mParameter.push_back(parameter{i});
    parameter{i}->units = units::unit_from_string("{pn['dimension']}");'''
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
EXPORTED RCPLIB::RCP<Resistance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        '''
        resistanceFriends.append(f'friend RCPLIB::RCP<Resistance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a Resistance instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<Resistance> create{name}(const RCPLIB::RCP< ElementImpl > &proxy){{
    auto ptr = __createResistance_do_not_call_outside_of_lib(proxy);
    ptr->numStates = 0;
    ptr->setDomain("{domain}");'''
        for i, param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    RCPLIB::RCP<Value> parameter{i} = RCPLIB::rcp(new Value("{param}", "{pn['value']}"));
    ptr->mParameter.push_back(parameter{i});
    parameter{i}->units = units::unit_from_string("{pn['dimension']}");'''
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
EXPORTED RCPLIB::RCP<PotentialSource> create{name}(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        '''
        PotentialSourceFriends.append(f'friend RCPLIB::RCP<PotentialSource> create{name}(const RCPLIB::RCP< ElementImpl > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a PotentialSource instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<PotentialSource> create{name}(const RCPLIB::RCP< ElementImpl > &proxy){{
    auto ptr = createPotentialSource(proxy);
    ptr->setDomain("{domain}");'''
        for i, param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    ptr->mParameter[{i}]->name = "{param}";
    ptr->mParameter[{i}]->value = "{pn['value']}";
    ptr->mParameter[{i}]->units = units::unit_from_string("{pn['dimension']}");'''
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
EXPORTED RCPLIB::RCP<FlowSource> create{name}(const RCPLIB::RCP< ElementImpl > &proxy=RCPLIB::null);

        '''
        flowSourceFriends.append(f'friend RCPLIB::RCP<FlowSource> create{name}(const RCPLIB::RCP< ElementImpl > &proxy);')
        hfile += decl
        decl = f'''
/**
* @brief Create a FlowSource instance, {description}
* Physical domain {domain}
* @return Reference counted pointer to the element 
*/
RCPLIB::RCP<FlowSource> create{name}(const RCPLIB::RCP< ElementImpl > &proxy){{
    auto ptr = createFlowSource(proxy);
    ptr->setDomain("{domain}");'''
        for i, param in enumerate(parameters.keys()):
            pn = parameters[param]
            isUC = "false"
            decl += f'''
    ptr->mParameter[{i}]->name = "{param}";
    ptr->mParameter[{i}]->value = "{pn['value']}";
    ptr->mParameter[{i}]->units = units::unit_from_string("{pn['dimension']}");'''
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
    

if __name__ == '__main__':
    spath = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(spath,'constitutiveRelations.json'),'r') as jfile:
        constitutiverelations = json.load(jfile)

        for domain,desc in constitutiverelations.items():
            generateStatementsForDomain(desc, domain)
        with open(os.path.join(spath,'../factorymethods.h'), 'w') as hf:
            print(hfile,file=hf)
        with open(os.path.join(spath,'../factorymethods.cpp'), 'w') as hc:
            print(cppfile,file=hc)

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

            if len(PotentialSourceFriends)>0:
                res = "#define DEFINE_FRIENDS_OF_PotentialSource \\"
                if len(PotentialSourceFriends) > 1:
                    for cp in PotentialSourceFriends[:-1]:
                        res += '\n\t'+cp+'\t\\'
                res += '\n\t'+PotentialSourceFriends[-1]
                stmts.append(res)
            else:
                stmts.append("#define DEFINE_FRIENDS_OF_PotentialSource")

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
            stmts.append("\n\nEXPORTED inline nlohmann::json getSupportedPhysicalDomainsAndFactoryMethods() {")
            stmts.append("\tnlohmann::json result;")
            dlist = '"Biochemical",'

            mlist = []
            for d,mts in createMethods.items():
                dlist += f'"{d}",'
                stmts.append(f"\t\n\tnlohmann::json mts{d};")
                for md,mn in mts.items():
                    stmts.append(f'\tmts{d}["{md}"]="{mn}";')
                    mlist.append(mn)
                stmts.append(f'\tresult["{d}"] = mts{d};\n\t')

            stmts.append(f'\tresult["Domains"] = nlohmann::json::array({{{dlist[:-1]}}});')
            stmts.append("\treturn result;")
            stmts.append("}")
            stmts.append('\n\n//Reversing strings to optimise string comparisions as most methods start with create')
            stmts.append("#define CALL_FACTORY_METHODS_BY_NAME \\")
            stmts.append("\tstd::string rMethodName(methodName.rbegin(),methodName.rend()); \\")
            stmts.append("\tRCPLIB::RCP<ElementImpl> data = RCPLIB::null;  \\")
            stmts.append(f'\tif(rMethodName=="{mlist[0][::-1]}") {{return {mlist[0]}(data);}} \\')
            for m in mlist[1:]:
                stmts.append(f'\telse if(rMethodName=="{m[::-1]}") {{return {m}(data);}} \\')

            stmts.append("\n#endif") #End header
            print('\n'.join(stmts),file=frd)