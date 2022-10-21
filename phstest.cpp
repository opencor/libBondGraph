#include "src/bondgraph.hpp"

using namespace BG;

RCPLIB::RCP<BondGraphInterface> rlc() {
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
  auto lR2 = createResistor();
  ioBondGraph->addComponent(lR2);

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
  ioBondGraph->connect(lJ0_1, lR2);
  ioBondGraph->connect(lSf, lJ0_1);

  std::cout << " PHS " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> reaction() {
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto A = createConcentration();
  A->setName("A");
  ioBondGraph->addComponent(A);

  auto B = createConcentration();
  B->setName("B");
  ioBondGraph->addComponent(B);

  // Create the reaction
  auto Re = createReaction();
  Re->setName("Re");
  ioBondGraph->addComponent(Re);

  // Create the junctions
  auto Y_A = createOneJunction();
  Y_A->setName("Y_A_1");
  ioBondGraph->addComponent(Y_A);

  auto Y_B = createOneJunction();
  Y_B->setName("Y_B_1");
  ioBondGraph->addComponent(Y_B);

  // Create the bonds
  ioBondGraph->connect(A, Y_A);
  ioBondGraph->connect(B, Y_B);
  ioBondGraph->connectInverting(Re, 0, Y_A);
  ioBondGraph->connectInverting(Re, 1, Y_B);

  std::cout << " Reaction " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> eReaction() {
  // Model reaction using electrical elements
  // Reaction is modelled as 1 junction with a resistor
  auto ioBondGraph = createBondGraph();

  // Create the storage
  auto lC1 = createCapacitor();
  lC1->setName("A");
  ioBondGraph->addComponent(lC1);

  auto lC2 = createCapacitor();
  lC2->setName("B");
  ioBondGraph->addComponent(lC2);

  // Create the resistor
  auto lR = createResistor();
  lR->setName("R");
  ioBondGraph->addComponent(lR);

  // Create the junctions
  auto lJ1_A = createOneJunction();
  auto lJ1_B = createOneJunction();
  auto lJ1_Re = createOneJunction();

  ioBondGraph->addComponent(lJ1_A);
  ioBondGraph->addComponent(lJ1_B);
  ioBondGraph->addComponent(lJ1_Re);

  // Create the bonds
  ioBondGraph->connect(lJ1_Re, lR);
  ioBondGraph->connect(lJ1_Re, lJ1_A);
  ioBondGraph->connect(lJ1_Re, lJ1_B);
  ioBondGraph->connect(lJ1_A, lC1);
  ioBondGraph->connect(lJ1_B, lC2);

  std::cout << " eReaction " << std::endl;
  ioBondGraph->computePortHamiltonian();

  return ioBondGraph;
}

int main(int argc, char *argv[]) {
  // rlc();
  reaction();
  // eReaction();
  return 0;
}