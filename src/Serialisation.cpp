#include "Exceptions.h"
#include "Serialisation.h"
#include "cellmlunitsmap.h"
#include "thirdparty/tinyxml2.h"
#include <sstream>
#include <symengine/eval.h>
#include <symengine/matrix.h>
#include <symengine/parser.h>
#include <symengine/parser/parser.h>
#include <symengine/solve.h>
#include <symengine/visitor.h>
#include <tuple>

namespace BG {
std::string getMathML(const SymEngine::RCP<const SymEngine::Basic> &expr)
{
    std::string expression = SymEngine::mathml(*expr);
    const std::string from = "type=\"integer\"";
    const std::string to = "cellml:units=\"dimensionless\"";

    auto toCellMLUnit = [&expression, from, to]() {
        size_t start_pos = 0;
        while ((start_pos = expression.find(from, start_pos)) != std::string::npos) {
            expression.replace(start_pos, from.length(), to);
            start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
        }
        return expression;
    };

    toCellMLUnit();

    return expression;
}

std::vector<std::string> getMathML(SymEngine::vec_basic &expressions)
{
    std::vector<std::string> mexpr;
    std::ostringstream ss;
    for (auto expr : expressions) {
        ss.str("");
        ss.clear();
        ss << "<apply><eq/>" << getMathML(expr) << "<cn cellml:units=\"dimensionless\">0.0</cn></apply>";
        std::string expression = ss.str();
        mexpr.push_back(expression);
    }
    return mexpr;
}

std::vector<std::string> getMathML(SymEngine::map_basic_basic &expressions)
{
    std::vector<std::string> mexpr;
    std::ostringstream ss;
    for (auto expr : expressions) {
        ss.str("");
        ss.clear();
        ss << *expr.first;

        std::string term = ss.str();
        ss.str("");
        ss.clear();
        ss << "<apply><eq/>" << getMathML(expr.first) << getMathML(expr.second) << "</apply>";
        mexpr.push_back(ss.str());
    }
    return mexpr;
}

/**
 * @brief Generate the cellml code for import and linking external component
 * 
 * @param vname Name of the variable that needs to be linked
 * @param prop Properties associated with that variable including the JSON description of the target file
 * import statement, connection statement, cellml statement (if import statement, this will be ""), variable name, bool (true if dot_vname is to be mapped)
 * @return std::tuple<std::string, std::string, std::string, std::string> 
 */
std::tuple<std::string, std::string, std::string, std::string, bool> generateLinkages(std::string vname, std::tuple<std::string, std::string, char> &prop)
{
    std::string value = std::get<1>(prop);
    try {
        nlohmann::json j = nlohmann::json::parse(value);
        if (std::get<2>(prop) == 'c') { //Control variable
            //Has vname and dot_vname
                std::string compartment = j["compartment"];
                std::string targetVar = vname;
                if (j.contains("target")) {
                    targetVar = j["target"];
                }
                std::string referenceName = vname;
                if (j.contains("importName")) {
                    referenceName = j["importName"];
                }
                std::string filename = j["filename"];
                std::string import = "<import xlink:href = \"" + filename + "\" >";
                import = import + "<component component_ref = \"" + compartment + "\" name = \"" + referenceName + "_IMP\" /> </import>";
                //std::string connection = "<connection>\n<map_components component_1 = \"" + vname + "_IMP\" component_2 = \"main\" />";
                //connection = connection + "<map_variables variable_1 = \""+targetVar+"\" variable_2 = \""+vname+"\" />";
                std::string connection = "<map_variables variable_1 = \"" + targetVar + "\" variable_2 = \"" + vname + "\" />";
                bool derivativeMapped = false;
                
                if (j.contains("mapvariables")) {
                    std::vector<std::string> vars = j["mapvariables"];
                    for (auto c : vars) {
                        connection = connection + "\n<map_variables variable_1 = \"" + c + "\" variable_2 = \"" + c + "\" />";
                    }
                }
                if(j.contains("derivative")){
                    derivativeMapped = true;
                    std::string dn = j["derivative"];
                    connection = connection + "\n<map_variables variable_1 = \"" + dn + "\" variable_2 = \"dot_" + vname + "\" />";
                }
                //connection = connection + "</connection>";
                return std::make_tuple(filename + ":" + import, compartment + ":" + referenceName + ":" + connection, "", vname,derivativeMapped);
        } else { //User defined
            std::string compartment = j["compartment"];
            std::map<std::string, std::string> targetVars = j["target"];
            std::string targetVar = targetVars[vname];
            std::string referenceName = j["importName"];
            std::string filename = j["filename"];
            std::string import = "<import xlink:href = \"" + filename + "\" >";
            import = import + "<component component_ref = \"" + compartment + "\" name = \"" + referenceName + "_IMP\" /> </import>";
            //std::string connection = "<connection>\n<map_components component_1 = \"" + vname + "_IMP\" component_2 = \"main\" />";
            //connection = connection + "<map_variables variable_1 = \"" + targetVar + "\" variable_2 = \"" + vname + "\" />";
            std::string connection = "<map_variables variable_1 = \"" + targetVar + "\" variable_2 = \"" + vname + "\" />";
            if (j.contains("mapvariables")) {
                std::vector<std::string> vars = j["mapvariables"];
                for (auto c : vars) {
                    connection = connection + "\n<map_variables variable_1 = \"" + c + "\" variable_2 = \"" + c + "\" />";
                }
            }
            bool derivativeMapped = false;
            if(j.contains("derivative")){
                derivativeMapped = true;
                std::map<std::string,std::string> dn = j["derivative"];
                for (auto c : dn) {                
                    connection = connection + "\n<map_variables variable_1 = \"" + c.first + "\" variable_2 = \"" + c.second + "\" />";
                }
            }            
            //connection = connection + "</connection>";
            return std::make_tuple(filename + ":" + import, compartment + ":" + referenceName + ":" + connection, "", vname,derivativeMapped);
        }
    } catch (std::exception &e) {
        throw BGException("Failed to generate cellml linkages for " + vname + " with value " + value + " (A valid JSON is expected)");
    }
}

std::map<std::string, std::string>
getCellML(std::string modelName, const BondGraph &host, std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>>> &bgequations)
{
    bool solvable;
    SymEngine::map_basic_basic equations;
    SymEngine::map_basic_basic bondEquations;
    SymEngine::vec_basic constraints;
    std::unordered_map<std::string, std::tuple<std::string, std::string, char>> dimensions;
    std::map<std::string, std::string> files;
    std::vector<std::string> connections;
    std::vector<std::string> imports;
    std::vector<std::string> blocks;

    std::ostringstream cellML;
    std::tie(solvable, equations, bondEquations, constraints, dimensions) = bgequations;
    //Find equations whose lhs has dot_
    std::unordered_map<std::string, std::string> derivatives;
    for (auto eq : equations) {
        cellML.str("");
        cellML.clear();
        cellML << *eq.first;
        std::string variable = cellML.str();
        if (variable.rfind("dot_", 0) == 0) {
            derivatives[variable.substr(4)] = variable;
        }
    }

    cellML.str("");
    cellML.clear();
    std::ostringstream mapParameterFile;
    std::ostringstream connect;

    mapParameterFile << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    mapParameterFile << "<model name=\"" << modelName << "_parameters\" xmlns=\"http://www.cellml.org/cellml/1.1#\" xmlns:cellml=\"http://www.cellml.org/cellml/1.1#\">" << std::endl;
    CellMLUnits unitMap;
    std::map<std::string, bool> mapUnits;
    std::map<std::string, std::string> defineUnits;
    std::map<std::string, std::string> definedNames;
    //Create variables and their annotations
    std::map<std::string, std::tuple<std::string, std::string, char>> orderedDim(dimensions.begin(), dimensions.end());
    for (auto var : orderedDim) {
        std::string dim = std::get<0>(var.second);
        std::string unitName = dim;
        try {
            unitName = unitMap.getUnitName(dim);
            definedNames[var.first] = unitName;
            mapUnits[unitName] = true;
        } catch (std::exception &x) {
            auto res = unitMap.getCellMLDef(unitName);
            defineUnits[unitName] = std::get<1>(res);
            definedNames[var.first] = std::get<0>(res);
        }
    }
    //Define cellml units
    //Map known units
    {
        tinyxml2::XMLDocument doc;
        tinyxml2::XMLPrinter printer;
        std::string units = unitMap.getPredinedunitsXML();
        if (doc.Parse(units.c_str()) == tinyxml2::XML_SUCCESS) {
            doc.Print(&printer);
            files["Units.cellml"] = printer.CStr();
        } else {
            files["Units.cellml"] = units;
            logWarn("Failed to parse standard units definition file!!");
        }
    }
    connect.str("");
    connect.clear();
    connect << "<import xlink:href = \"Units.cellml\">" << std::endl;
    mapParameterFile << "<import xlink:href = \"Units.cellml\">" << std::endl;
    for (auto mu : mapUnits) {
        connect << "<units name = \"" << mu.first << "\" units_ref=\"" << mu.first << "\" />" << std::endl;
        mapParameterFile << "<units name = \"" << mu.first << "\" units_ref=\"" << mu.first << "\" />" << std::endl;
    }
    connect << "</import>" << std::endl;
    mapParameterFile << "</import>" << std::endl;
    {
        std::string reS = connect.str();
        connections.push_back(reS); //Units
    }
    connect.str("");
    connect.clear();
    //Setup units in both files
    for (auto un : defineUnits) {
        connect << un.second << std::endl;
        mapParameterFile << un.second << std::endl;
    }
    {
        std::string reS = connect.str();
        blocks.push_back(reS); //Units
    }
    connect.str("");
    connect.clear();
    //Define link to parameters file
    connect << "<import xlink:href = \""
            << modelName << "_parameters.cellml\" >" << std::endl;
    connect << "<component component_ref = \"parameters\" name = \"parameters\" />" << std::endl;
    connect << "</import>" << std::endl;
    {
        std::string reS = connect.str();
        imports.push_back(reS); //Units
    }
    connect.str("");
    connect.clear();

    //Define all the parameters and map them to cellML
    connect << "<connection>\n<map_components component_1 = \"parameters\" component_2 = \"main\" />" << std::endl;

    mapParameterFile << "<component name=\"parameters\">" << std::endl;
    for (auto var : orderedDim) {
        if (std::get<2>(var.second) == 'p') {
            std::string dim = std::get<0>(var.second);
            std::string unitName = definedNames[var.first];
            connect << "<map_variables variable_1 = \"" << var.first << "\" variable_2 = \"" << var.first << "\" />" << std::endl;
            mapParameterFile << "<variable initial_value = \"" << std::get<1>(var.second) << "\" name = \"" << var.first << "\" public_interface = \"out\" units = \"" << unitName << "\" />" << std::endl;
        }
    }
    connect << "</connection>" << std::endl;
    mapParameterFile << "</component>\n</model>" << std::endl;
    {
        std::string reS = connect.str();
        connections.push_back(reS); //Units
    }
    {
        tinyxml2::XMLDocument doc;
        tinyxml2::XMLPrinter printer;
        std::string params = mapParameterFile.str();
        if (doc.Parse(params.c_str()) == tinyxml2::XML_SUCCESS) {
            doc.Print(&printer);
            files[modelName + "_parameters.cellml"] = printer.CStr();
        } else {
            files[modelName + "_parameters.cellml"] = params;
            logWarn("Failed to parse generated parameters definition file!!");
        }
    }
    std::unordered_map<std::string, bool> computedVariables;
    std::unordered_map<std::string, std::string> mappedVariables;
    for (auto eq : equations) {
        std::string eqs = eq.first->__str__();
        computedVariables[eqs] = true;
    }
    for (auto eq : bondEquations) {
        std::string eqs = eq.first->__str__();
        computedVariables[eqs] = true;
    }
    //All variables that appear in a constraint are computed i.e. no initial value should be provided
    for (auto c : constraints) {
        SymEngine::set_basic atoms = SymEngine::atoms<SymEngine::FunctionSymbol, SymEngine::Symbol>(*c);
        for (auto eq : atoms) {
            std::string eqs = eq->__str__();
            computedVariables[eqs] = true;
        }
    }
    //Find all variables that need to be mapped
    //Control variables, userdefined variables - values must be JSON
    std::map<std::string, std::string> importMap;
    std::map<std::string, std::vector<std::tuple<std::string, std::string>>> connectionMap;
    for (auto var : orderedDim) {
        std::string value = std::get<1>(var.second);
        //Check if json
        if (value.find('{') != std::string::npos) {
            //Generate import, connection and file for value and the name to import
            auto g = generateLinkages(var.first, var.second);
            std::string imp = std::get<0>(g);
            auto ss = imp.find(':');
            auto importFile = imp.substr(0, ss);
            importMap[importFile] = imp.substr(ss + 1);
            imp = std::get<1>(g);
            ss = imp.find(':');
            auto compartment = imp.substr(0, ss);
            auto nxbit = imp.substr(ss+1);
            auto si = nxbit.find(':');
            auto variable = nxbit.substr(0,si);

            if (connectionMap.find(compartment) != connectionMap.end()) {
                connectionMap[compartment].push_back(std::make_tuple(variable, nxbit.substr(si + 1)));
            } else {
                connectionMap[compartment] = {std::make_tuple(variable, nxbit.substr(si + 1))};
            }
            if (std::get<2>(g) != "")
                files[var.first + ".cellml"] = std::get<2>(g);
            mappedVariables[var.first] = std::get<3>(g);
            if (std::get<2>(var.second) == 'c' && std::get<4>(g)) { //Create a map only if a mapping is available
                mappedVariables["dot_" + var.first] = "dot_" + std::get<3>(g);
            }
        }
    }
    for (auto c : importMap) {
        imports.push_back(c.second);
    }
    for (auto c : connectionMap) {
        std::map<std::string, std::vector<std::string>> vmap;
        for (auto v : c.second) {
            if (vmap.find(std::get<0>(v)) != vmap.end()) {
                vmap[std::get<0>(v)].push_back(std::get<1>(v));
            } else {
                vmap[std::get<0>(v)] = {std::get<1>(v)};
            }
        }
        std::string comp = c.first;
        for (auto v : vmap) {
            std::string connection = "<connection>\n<map_components component_1 = \"" + v.first + "_IMP\" component_2 = \"main\" />";
            for (auto x : v.second) {
                connection = connection + x;
            }
            connection = connection + "</connection>";
            connections.push_back(connection);
        }
    }
    //Start generation
    cellML.str("");
    cellML.clear();
    cellML << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    cellML << "<model name=\"" << modelName << "\" xmlns=\"http://www.cellml.org/cellml/1.1#\" xmlns:cellml=\"http://www.cellml.org/cellml/1.1#\">" << std::endl;
    //Define units
    cellML << connections[0] << std::endl;
    connections.erase(connections.begin()); //Remove units import
    cellML << blocks[0] << std::endl;
    blocks.erase(blocks.begin()); //Remove units
    //Define imports
    for (auto c : imports) {
        cellML << c << std::endl;
    }
    //Define connections
    for (auto c : connections) {
        cellML << c << std::endl;
    }
    //Define other blocks
    for (auto c : blocks) {
        cellML << c << std::endl;
    }
    //Define the component
    cellML << "<component name=\"main\">" << std::endl;

    if (!solvable) {
        cellML << "<!-- The bondgraph library has identified issues with this model. The cellml description may not build or give correct results, view the library's log (if enabled) for more details -->" << std::endl;
    }
    //Define time variable
    cellML << "<variable name=\"t\" units=\"second\"/>" << std::endl;

    for (auto var : orderedDim) {
        std::string dim = std::get<0>(var.second);
        std::string unitName = definedNames[var.first];
        std::string initVal = "initial_value = \"" + std::get<1>(var.second) + "\"";
        if (computedVariables[var.first]) {
            initVal = "";
        }
        switch (std::get<2>(var.second)) {
        case 'p': {
            cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName << "\" public_interface = \"in\" />";
            break;
        }
        case 'c': {
            //TODO: These need to be mapped to a component
            if(mappedVariables.find(var.first)==mappedVariables.end()){
                cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName << "\" " << initVal << "/>";
            }else{
                cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName << "\" public_interface = \"in\" />";
            }
            break;
        }
        case 's': {
            if (var.first.size() > 4 && var.first.substr(0, 4) == "dot_") { //Derivatives are computed
                cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName << "\" public_interface = \"out\" />";
            } else {
                cellML << "< variable  name=\"" << var.first << "\" units=\"" << unitName << "\" public_interface = \"out\" " << initVal << "/>";
            }
            break;
        }
        default: { //Variables
            cellML << "< variable name=\"" << var.first << "\" units=\"" << unitName << "\" public_interface = \"out\" />";
            break;
        }
        }
    }
    cellML << "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">" << std::endl;
    //Create Bond variable equations
    std::vector<std::string> bondEqMathML = getMathML(bondEquations);
    cellML << "<!--Port values -->" << std::endl;
    for (auto expr : bondEqMathML) {
        cellML << expr << std::endl;
    }
    //Create Constraint equations
    if (constraints.size() > 0) {
        cellML << "<!--Constraint equations -->" << std::endl;
        std::vector<std::string> consEqMathML = getMathML(constraints);
        for (auto expr : consEqMathML) {
            cellML << expr << std::endl;
        }
    }
    std::vector<std::string> stateEqMathML = getMathML(equations);
    cellML << "<!--State equations -->" << std::endl;
    for (auto expr : stateEqMathML) {
        cellML << expr << std::endl;
    }
    for (auto expr : derivatives) {
        cellML << "<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>" << expr.first << "</ci></apply><ci>" << expr.second << "</ci></apply> " << std::endl;
    }
    cellML << "</math>" << std::endl;
    cellML << "</component>" << std::endl;
    cellML << "</model>";
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLPrinter printer;

    if (doc.Parse(cellML.str().c_str()) == tinyxml2::XML_SUCCESS) {
        doc.Print(&printer);
        files[modelName + ".cellml"] = printer.CStr();
        return files;
    } else {
        logDebug(cellML.str());
        throw BGException("Error creating CellML");
    }
}
} // namespace BG