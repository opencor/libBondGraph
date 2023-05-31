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
#include <gtest/gtest.h>

#include "bondgraph.hpp"
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <cstdio>

using namespace BG;
// Demonstrate some basic assertions.

RCPLIB::RCP<BondGraphInterface> rlc() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lI = createInductor();
  ioBondGraph->addComponent(lI);

  auto lC1 = createCapacitor();
  lC1->setParameter("C","1.0","mega Farad");
  ioBondGraph->addComponent(lC1);

  // Create the resistor
  auto lR = createResistor();
  lR->setParameter("R","0.5","micro Ohm");
  ioBondGraph->addComponent(lR);

  auto lE = createConstantVoltageSource();
  ioBondGraph->addComponent(lE);

  // Create the junctions
  // auto lJ0_1 = createZeroJunction();
  auto lJ0_1 = createOneJunction();
  ioBondGraph->addComponent(lJ0_1);

  // Create the bonds
  ioBondGraph->connect(lR, lJ0_1);
  ioBondGraph->connect(lI, lJ0_1);
  ioBondGraph->connect(lC1, lJ0_1);
  ioBondGraph->connect(lE, lJ0_1);

  auto res = ioBondGraph->computeStateEquation();
  std::cout << " RLC E " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> reaction() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto A = createConcentration();
  ioBondGraph->addComponent(A);

  auto B = createConcentration();
  ioBondGraph->addComponent(B);

  // Create the reaction
  auto Re = createReaction();
  ioBondGraph->addComponent(Re);

  // Create the junctions
  auto Y_A = createOneJunction();
  ioBondGraph->addComponent(Y_A);

  auto Y_B = createOneJunction();
  ioBondGraph->addComponent(Y_B);

  // Create the bonds
  ioBondGraph->connect(A, Y_A);
  ioBondGraph->connect(B, Y_B);
  ioBondGraph->connectInverting(Re, 0, Y_A);
  ioBondGraph->connectInverting(Re, 1, Y_B);

  auto res = ioBondGraph->computeStateEquation();
  std::cout << " Reaction " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> phstest() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createCapacitor();
  ioBondGraph->addComponent(lC1);

  auto lC2 = createCapacitor();
  ioBondGraph->addComponent(lC2);

  // Create the Flow source
  auto lSf = createConstantFlowSource();
  ioBondGraph->addComponent(lSf);

  // Create the resistor
  auto lR = createResistor();
  ioBondGraph->addComponent(lR);

  // Create the Transformer
  auto lTf = createTransformer();
  ioBondGraph->addComponent(lTf);

  // Create the junctions
  auto lJ0_1 = createZeroJunction();
  auto lJ1_1 = createOneJunction();
  ioBondGraph->addComponent(lJ0_1);
  ioBondGraph->addComponent(lJ1_1);

  // Create the bonds
  ioBondGraph->connect(lJ1_1, lR);
  ioBondGraph->connect(lJ1_1, lC1);
  ioBondGraph->connect(lJ1_1, lTf);

  ioBondGraph->connect(lTf, 1, lJ0_1);
  ioBondGraph->connect(lJ0_1, lC2);
  ioBondGraph->connect(lSf, lJ0_1);

  //std::cout << " PHS " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> phstestRLC() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createCapacitor();
  ioBondGraph->addComponent(lC1);

  auto lL2 = createInductor();
  ioBondGraph->addComponent(lL2);

  // Create the Flow source
  auto lSfin = createConstantFlowSource();
  ioBondGraph->addComponent(lSfin);

  auto lSfout = createConstantFlowSource();
  ioBondGraph->addComponent(lSfout);  

  // Create the resistor
  auto lR = createResistor();
  ioBondGraph->addComponent(lR);

  // Create the junctions
  auto lJ0_1 = createZeroJunction();
  auto lJ0_3 = createZeroJunction();  
  auto lJ1_2 = createOneJunction();
  ioBondGraph->addComponent(lJ0_1);
  ioBondGraph->addComponent(lJ0_3);
  ioBondGraph->addComponent(lJ1_2);

  // Create the bonds
  ioBondGraph->connect(lJ0_1, lC1);
  ioBondGraph->connect(lJ0_3, lL2);
  ioBondGraph->connect(lJ1_2, lR);

  ioBondGraph->connect(lSfin,lJ0_1);
  ioBondGraph->connect(lJ0_3, lSfout);

  ioBondGraph->connect(lJ0_1,lJ1_2);
  ioBondGraph->connect(lJ1_2,lJ0_3);

  auto res = ioBondGraph->computeStateEquation();
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> phstestHydraulic() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createFluidCompliance();
  ioBondGraph->addComponent(lC1);

  auto lL2 = createFluidInertance();
  ioBondGraph->addComponent(lL2);

  // Create the Flow source
  auto lSfin = createConstantFluidFlowSource();
  ioBondGraph->addComponent(lSfin);

  auto lSfout = createConstantFluidFlowSource();
  ioBondGraph->addComponent(lSfout);  

  // Create the resistor
  auto lR = createViscousResistance();
  ioBondGraph->addComponent(lR);

  // Create the junctions
  auto lJ0_1 = createZeroJunction();
  auto lJ0_3 = createZeroJunction();  
  auto lJ1_2 = createOneJunction();
  ioBondGraph->addComponent(lJ0_1);
  ioBondGraph->addComponent(lJ0_3);
  ioBondGraph->addComponent(lJ1_2);

  // Create the bonds
  ioBondGraph->connect(lJ0_1, lC1);
  ioBondGraph->connect(lJ0_3, lL2);
  ioBondGraph->connect(lJ1_2, lR);

  ioBondGraph->connect(lSfin,lJ0_1);
  ioBondGraph->connect(lJ0_3, lSfout);

  ioBondGraph->connect(lJ0_1,lJ1_2);
  ioBondGraph->connect(lJ1_2,lJ0_3);

  auto res = ioBondGraph->computeStateEquation();
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

TEST(BondGraphSetup, PHS) {
  newWorkSpace();
  // Expect two strings not to be equal.
  auto ioBondGraph = phstestRLC();
  EXPECT_EQ(ioBondGraph.is_null(), false);
}

TEST(BondGraphSetup, Electrical) {
  newWorkSpace();
  // Expect two strings not to be equal.
  auto ioBondGraph = phstestRLC();

  auto res = ioBondGraph->computeStateEquation();
  // Solvable
  EXPECT_EQ(res.bondGraphValidity, true);
  auto equations = res.dof;
  std::map<std::string, std::string> result = {
      {"dot_Inductor", "Capacitor/C_0 + r_2*i_4 - r_2*Inductor/L_1"},
      {"dot_Capacitor", "i_3 + i_4 - Inductor/L_1"}};

  for (auto eq : equations) {
    std::string svar = symEngineExpressionToString(eq.first);
    std::string sres = symEngineExpressionToString(eq.second);

    std::string expected = result[svar];
    EXPECT_EQ(sres, expected);
  }
}

TEST(BondGraphSetup, Hydraulic) {
  newWorkSpace();
  // Expect two strings not to be equal.
  auto ioBondGraph = phstestHydraulic();

  auto res = ioBondGraph->computeStateEquation();
  // Solvable
  EXPECT_EQ(res.bondGraphValidity, true);
  auto equations = res.dof;
  std::map<std::string, std::string> result = {
      {"dot_FluidCompliance", "i_3 + i_4 - FluidInertance/L_1"},
      {"dot_FluidInertance", "FluidCompliance/C_0 + r_2*i_4 - r_2*FluidInertance/L_1"}};
  for (auto eq : equations) {
    std::string svar = symEngineExpressionToString(eq.first);
    std::string sres = symEngineExpressionToString(eq.second);
    std::string expected = result[svar];
    EXPECT_EQ(sres, expected);
  }
}

TEST(BondGraphSetup, Biochemical) {
  newWorkSpace();
  auto ioBondGraph = reaction();

  auto res = ioBondGraph->computeStateEquation();
  // Solvable
  EXPECT_EQ(res.bondGraphValidity, true);
  auto equations = res.dof;
  std::map<std::string, std::string> result = {
      {"dot_Ce", "-Ce*r_2*k_0 + Ce*r_2*k_1"}};

  for (auto eq : equations) {
    std::string svar = symEngineExpressionToString(eq.first);
    std::string sres = symEngineExpressionToString(eq.second);
    std::string expected = result[svar];
    EXPECT_EQ(sres, expected);
  }
}

TEST(BondGraphSetup, Composition) {
  newWorkSpace();
  auto bg = reaction();
  auto parent = createBondGraph();
  parent->addBondGraph(bg);
  auto clone = bg->clone();
  parent->addBondGraph(clone);
  parent->computeStateEquation();
  auto res = parent->computeStateEquation();
  // Solvable
  EXPECT_EQ(res.bondGraphValidity, true);
  auto equations = res.dof;
  std::map<std::string, std::string> result = {
      {"dot_Ce", "-Ce*r_5*k_2 + Ce*r_5*k_3"},
      {"dot_a_1", "-(-k*a_0*r_5 + k*a_1*r_5)"},
      {"dot_a_2", "-(k*a_2*r_5 - k*a_3*r_5)"},
      {"dot_a_3", "-(-k*a_2*r_5 + k*a_3*r_5)"}};

  for (auto eq : equations) {
    std::string svar = symEngineExpressionToString(eq.first);
    std::string sres = symEngineExpressionToString(eq.second);
    std::string expected = result[svar];
    EXPECT_EQ(sres, expected);
  }
}

TEST(BondGraphSetup, Serialisation) {
  newWorkSpace();
  auto bg = reaction();
  auto parent = createBondGraph();
  parent->addBondGraph(bg);
  auto clone = bg->clone();
  clone->computeStateEquation();
  parent->addBondGraph(clone);

  auto eqs = bg->computeStateEquation();
  auto files = getCellML("Reaction", bg, eqs);
  char  fname[256];
	std::tmpnam(fname);  
  for (const auto &itm : files) {
    //std::string fname = "D:/Temp/" + itm.first;
    std::ofstream ofs(fname, std::ofstream::out);
    ofs << itm.second;
    ofs.close();
  }

  nlohmann::json j = parent;

  newWorkSpace();

  RCPLIB::RCP<BondGraphInterface> bg2 = createBondGraph();
  j.get_to(bg2);

  nlohmann::json j2 = bg2;
  EXPECT_EQ((j == j2), true);
}

TEST(BondGraphSetup, CellML) {
  newWorkSpace();
  auto ioBondGraph = createBondGraph();

  // Create the resistor
  auto lR = createResistor();
  lR->setParameter("r", "1", "Ohm");
  ioBondGraph->addComponent(lR);

  auto lC1 = createCapacitor();
  lC1->setParameter("C", "1", "Farad");
  lC1->setParameter("q_0", "0", "coulomb");
  ioBondGraph->addComponent(lC1);

  // Create the junctions
  auto lJ0_1 = createOneJunction();
  ioBondGraph->addComponent(lJ0_1);

  // Create the bonds
  ioBondGraph->connect(lR, lJ0_1);
  ioBondGraph->connect(lC1, lJ0_1);

  auto eqs = ioBondGraph->computeStateEquation();
  auto files = getCellML("RLC", ioBondGraph, eqs);
  EXPECT_EQ(true, files.find("RLC.cellml")!=files.end());
  EXPECT_EQ(true, files.find("RLCParameters.cellml")!=files.end());
  EXPECT_EQ(true, files.find("Units.cellml")!=files.end());
}

TEST(BondGraphSetup, SupportedDomains) {
  newWorkSpace();
  auto sd = getSupportedPhysicalDomainsAndFactoryMethods();
  EXPECT_EQ(true, sd.contains("ElementDefinitions"));
}
