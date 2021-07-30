#ifndef SERIALIZATION_H__
#define SERIALIZATION_H__
#include <vector>
#include <string>
#include <tuple>
#include "bondgraph.h"
#include <symengine/printers.h>

namespace BG{

    std::string getMathML(const SymEngine::RCP<const SymEngine::Basic>& expr);

    std::vector<std::string> getMathML(SymEngine::vec_basic& expressions);

    std::vector<std::string> getMathML(SymEngine::map_basic_basic& expressions);

    std::map<std::string, std::string> getCellML(std::string modelName, const BondGraph &host, std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>> > &bgequations);
}

#endif