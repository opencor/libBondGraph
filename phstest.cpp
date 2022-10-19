#include "src/bondgraph.hpp"

using namespace BG;

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

int main(int argc, char *argv[]) {
  phstest();
  return 0;
}