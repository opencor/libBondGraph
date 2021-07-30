#include <gtest/gtest.h>

#include "Serialisation.h"
#include "bondgraph.h"
#include "componentregistry.h"
#include <sstream>
#include <string>
#include <symengine/basic.h>
#include <tuple>
#include <vector>

using namespace BG;
// Demonstrate some basic assertions.

RCPLIB::RCP<BondGraph> rlc()
{
    auto ioBondGraph = createBondGraph();

    //Create the storage
    auto lI = createInductor();
    ioBondGraph->addComponent(lI);

    auto lC1 = createCapacitor();
    ioBondGraph->addComponent(lC1);

    //Create the resistor
    auto lR = createResistor();
    ioBondGraph->addComponent(lR);

    auto lE = createConstantVoltageSource();
    ioBondGraph->addComponent(lE);

    //Create the junctions
    auto lJ0_1 = createZeroJunction();
    ioBondGraph->addComponent(lJ0_1);

    //Create the bonds
    ioBondGraph->connect(lR, lJ0_1);
    ioBondGraph->connect(lI, lJ0_1);
    ioBondGraph->connect(lC1, lJ0_1);
    ioBondGraph->connect(lE, lJ0_1);

    return ioBondGraph;
}

RCPLIB::RCP<BondGraph> reaction()
{
    auto ioBondGraph = createBondGraph();

    //Create the storage
    auto A = createConcentration();
    ioBondGraph->addComponent(A);

    auto B = createConcentration();
    ioBondGraph->addComponent(B);

    //Create the reaction
    auto Re = createReaction();
    ioBondGraph->addComponent(Re);

    //Create the junctions
    auto Y_A = createOneJunction();
    ioBondGraph->addComponent(Y_A);

    auto Y_B = createOneJunction();
    ioBondGraph->addComponent(Y_B);

    //Create the bonds
    ioBondGraph->connect(A, Y_A);
    ioBondGraph->connect(B, Y_B);
    ioBondGraph->connectInverting(Re, 0, Y_A);
    ioBondGraph->connectInverting(Re, 1, Y_B);
    return ioBondGraph;
}

TEST(BondGraphSetup, Electrical)
{
    newWorkSpace();
    // Expect two strings not to be equal.
    auto ioBondGraph = rlc();

    std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>>> res = ioBondGraph->computeStateEquation();
    //Solvable
    EXPECT_EQ(std::get<0>(res), true);
    auto equations = std::get<1>(res);
    std::map<std::string, std::string> result = {
        {"q_1","u_3*C_1"},
        {"dot_p_0","u_3"},
        {"dot_q_1","C_1*dot_u_3"}
    };
    for (auto eq : equations) {
        std::string svar = eq.first->__str__();
        std::string sres = eq.second->__str__();
        std::string expected = result[svar];
        EXPECT_EQ(sres, expected);
    }
}

TEST(BondGraphSetup, Biochemical)
{
    newWorkSpace();
    auto ioBondGraph = reaction();

    std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>>> res = ioBondGraph->computeStateEquation();
    //Solvable
    EXPECT_EQ(std::get<0>(res), true);
    auto equations = std::get<1>(res);
    std::map<std::string, std::string> result = {
        {"dot_a_0" , "-(k*a_0*r_2 - k*a_1*r_2)"},
        {"dot_a_1" , "-(-k*a_0*r_2 + k*a_1*r_2)"}};
    for (auto eq : equations) {
        std::string svar = eq.first->__str__();
        std::string sres = eq.second->__str__();
        std::string expected = result[svar];
        //std::cout<<"{\""<<svar<<"\",\""<<sres<<"\"}"<<std::endl;
        EXPECT_EQ(sres, expected);
    }
}

TEST(BondGraphSetup, Composition)
{
    newWorkSpace();
    auto bg = reaction();
    auto parent = createBondGraph();
    parent->addBondGraph(bg);
    auto clone = bg->clone();
    parent->addBondGraph(clone);
    parent->computeStateEquation();
    std::tuple<bool, SymEngine::map_basic_basic, SymEngine::map_basic_basic, SymEngine::vec_basic, std::unordered_map<std::string, std::tuple<std::string, std::string, char>>> res = parent->computeStateEquation();
    //Solvable
    EXPECT_EQ(std::get<0>(res), true);
    auto equations = std::get<1>(res);
    std::map<std::string, std::string> result = {
        {"dot_a_0","-(k*a_0*r_5 - k*a_1*r_5)"},
        {"dot_a_1","-(-k*a_0*r_5 + k*a_1*r_5)"},
        {"dot_a_2","-(k*a_2*r_5 - k*a_3*r_5)"},
        {"dot_a_3","-(-k*a_2*r_5 + k*a_3*r_5)"}
    };

    for (auto eq : equations) {
        std::string svar = eq.first->__str__();
        std::string sres = eq.second->__str__();
        std::string expected = result[svar];
        EXPECT_EQ(sres, expected);
    }
}

TEST(BondGraphSetup, Serialisation)
{
    newWorkSpace();
    auto bg = reaction();
    auto parent = createBondGraph();
    parent->addBondGraph(bg);
    auto clone = bg->clone();
    clone->computeStateEquation();
    parent->addBondGraph(clone);

    nlohmann::json j = parent;

    newWorkSpace();

    RCPLIB::RCP<BondGraph> bg2 = createBondGraph();
    j.get_to(bg2);

    nlohmann::json j2 = bg2;
    EXPECT_EQ((j == j2),true);
}

TEST(BondGraphSetup,CellML){
    newWorkSpace();
    auto ioBondGraph = createBondGraph();

    //Create the storage
    auto lI = createInductor();
    ioBondGraph->addComponent(lI);

    auto lC1 = createCapacitor();
    ioBondGraph->addComponent(lC1);

    //Create the resistor
    auto lR = createResistor();
    ioBondGraph->addComponent(lR);

    auto lE = createConstantVoltageSource();
    lE->setValue(std::string("u"), std::string("{\"type\": \"file\",\"compartment\": \"cvolt\",\"filename\": \"u1.cellml\",\"mapvariables\": [\"t\"],\"derivative\":\"du\"}"));
    ioBondGraph->addComponent(lE);

    //Create the junctions
    auto lJ0_1 = createZeroJunction();
    ioBondGraph->addComponent(lJ0_1);

    //Create the bonds
    ioBondGraph->connect(lR, lJ0_1);
    ioBondGraph->connect(lI, lJ0_1);
    ioBondGraph->connect(lC1, lJ0_1);
    ioBondGraph->connect(lE, lJ0_1);

    auto eqs = ioBondGraph->computeStateEquation();
    auto files = getCellML("RLC", *ioBondGraph, eqs);
    std::map<std::string, int> expected = {
        {"RLC.cellml", 11363},
        {"RLC_parameters.cellml", 1601}};
    for (auto c : expected) {
        int rr = files[c.first].size();
        EXPECT_EQ(rr,c.second);
    }
}