#include "bondgraph.h"
#include "Serialisation.h"
#include "componentregistry.h"
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace BG;

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

    /*
    SymEngine::vec_basic equations;
    SymEngine::vec_basic bondEquations;
    SymEngine::vec_basic constraints;

    std::tie(equations, bondEquations, constraints) = ioBondGraph->computeStateEquation();
    std::cout << equations << std::endl;
    std::cout << constraints << std::endl;
    */
    ioBondGraph->computeStateEquation();
    //std::cout<<ioBondGraph->serializeAsCellML("RLC")<<std::endl;
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

    //auto Y_A_proxy = createProxy(Y_A);
    //ioBondGraph->addComponent(Y_A_proxy);

    //Create the bonds
    ioBondGraph->connect(A, Y_A);
    ioBondGraph->connect(B, Y_B);
    ioBondGraph->connectInverting(Re, 0, Y_A);
    ioBondGraph->connectInverting(Re, 1, Y_B);

    /*
    SymEngine::vec_basic equations;
    SymEngine::vec_basic bondEquations;
    SymEngine::vec_basic constraints;

    std::tie(equations, bondEquations, constraints) = ioBondGraph->computeStateEquation();
    std::cout << equations << std::endl;
    std::cout << constraints << std::endl;
    */
    ioBondGraph->computeStateEquation();
    //std::cout<<ioBondGraph->serializeAsCellML("Reaction")<<std::endl;
    return ioBondGraph;
}

void checkAdditions()
{
    std::cerr << "Reaction " << std::endl;
    auto bg = reaction();
    bg->computeStateEquation();
    auto parent = createBondGraph();
    parent->addBondGraph(bg);
    parent->computeStateEquation();
    std::cerr << "Clone" << std::endl;
    auto clone = bg->clone();
    clone->computeStateEquation();
    std::cerr << "Clone add" << std::endl;
    parent->addBondGraph(clone);
    parent->computeStateEquation();
    {
        std::cout << "\n";
        std::cout << "BG components " << std::endl;
        auto bgComps = bg->getComponents();
        for (auto c : bgComps) {
            std::cout << std::get<0>(c) << "\t is proxy " << std::get<2>(c) << std::endl;
        }
        std::cout << "Parent components " << std::endl;
        auto parentComps = parent->getComponents();
        for (auto c : parentComps) {
            std::cout << std::get<0>(c) << "\t is proxy " << std::get<2>(c) << std::endl;
        }
        std::cout << "Name map " << std::endl;
        auto nm = ComponentRegistry::defaultRegistry().getNameMap();
        for (auto n : nm) {
            std::cout << n.first << "\t" << n.second << std::endl;
        }
    }
    std::cout << "\nRemoving bondgraph " << std::endl;
    parent->removeBondGraph(clone);
    parent->computeStateEquation();
    {
        std::cout << "\n";
        std::cout << "BG components " << std::endl;
        auto bgComps = bg->getComponents();
        for (auto c : bgComps) {
            std::cout << std::get<0>(c) << "\t is proxy " << std::get<2>(c) << std::endl;
        }
        std::cout << "Parent components " << std::endl;
        auto parentComps = parent->getComponents();
        for (auto c : parentComps) {
            std::cout << std::get<0>(c) << "\t is proxy " << std::get<2>(c) << std::endl;
        }
        std::cout << "Name map " << std::endl;
        auto nm = ComponentRegistry::defaultRegistry().getNameMap();
        for (auto n : nm) {
            std::cout << n.first << "\t" << n.second << std::endl;
        }
    }
    //rlc();
}


void portHamiltonian(){
    auto ioBondGraph = createBondGraph();
    auto port_hamiltonian = createPortHamiltonian("(w +G*x_0)*(x_1**2+x_2**2) / 2", {"x_0", "x_1", "x_2"});
    //Set SI units - by default they are set to dimensionless
    /*
    port_hamiltonian->setSIUnit("w","");
    port_hamiltonian->setSIUnit("G", "");
    port_hamiltonian->setSIUnit("x_0", "");
    port_hamiltonian->setSIUnit("x_1", "");
    port_hamiltonian->setSIUnit("x_2", "");
    */
   
    ioBondGraph->addComponent(port_hamiltonian);
    //Define the symplectic junction structure
    auto symplectic_gyrator = createGyrator();
    ioBondGraph->addComponent(symplectic_gyrator);
    auto em_field = createOneJunction();
    ioBondGraph->addComponent(em_field);

    ioBondGraph->connect(port_hamiltonian,1,em_field);
    ioBondGraph->connect(symplectic_gyrator, 1, em_field);
    ioBondGraph->connect(port_hamiltonian,2,symplectic_gyrator,0);
    
    //Construct the open part of the system
    auto dissipation = createResistor();
    ioBondGraph->addComponent(dissipation);

    auto photon_source = createConstantCurrentSource();
    ioBondGraph->addComponent(photon_source);
    
    ioBondGraph->connect(em_field, dissipation);
    ioBondGraph->connect(photon_source, em_field);
    ioBondGraph->computeStateEquation();
}

#include "../src/thirdparty/json.hpp"
#include "../src/ElementsImpl.h"

void testJSON(){
    auto bg = reaction();
    auto parent = createBondGraph();
    parent->addBondGraph(bg);
    auto clone = bg->clone();    
    //clone->computeStateEquation();
    parent->addBondGraph(clone);
    //parent->computeStateEquation();
    
    nlohmann::json j = parent;

    //std::cout << j << std::endl;
    newWorkSpace();

    RCPLIB::RCP<BondGraph> bg2 = createBondGraph();
    j.get_to(bg2);

    nlohmann::json j2 = bg2;
    std::cout<<" Compare "<< (j==j2) <<std::endl;
    std::cout<<"First\n\n"<<j<<"\nSecond\n\n"<<j2<<std::endl;
    //loadJson(j);
}

#include <units.hpp>

int main(int argc, char *argv[])
{
    //checkAdditions();
    //testJSON();
    //portHamiltonian();

    auto bg = rlc();
    auto eqs = bg->computeStateEquation();
    auto files = getCellML("RLC",*bg,eqs);
    std::cout<<files["RLC.cellml"]<<std::endl;

    /*
    rlc();
    newWorkSpace();
    reaction();
    newWorkSpace();
    
    newWorkSpace();
    portHamiltonian();

    std::cout << getSupportedPhysicalDomainsAndFactoryMethods() << std::endl;
    auto cap = createBondgraphElement("createCapacitor");
    std::cout<<" Type "<<cap->getType()<<"\t"<<cap->getDomain()<<std::endl;
    std::cerr
        << "Reaction " << std::endl;
    auto bg = reaction();
    //bg->computeStateEquation();
    std::cerr << "RLC " << std::endl;
    auto bgrlc = rlc();
    //bgrlc->computeStateEquation();
    std::cerr << "Parent " << std::endl;
    auto parent = createBondGraph();
    parent->addBondGraph(bg);
    parent->addBondGraph(bgrlc);
    //parent->computeStateEquation();
    */
    /*
    auto nm = ComponentRegistry::defaultRegistry().getNameMap();
    for(auto n : nm){
        std::cout<<n.first<<"\t"<<n.second<<std::endl;
    } 
    */
    /*
    auto &OJ2 = std::get<1>(ComponentRegistry::defaultRegistry().getComponentByName("1J_1"));
    auto &Se = std::get<1>(ComponentRegistry::defaultRegistry().getComponentByName("Se_1"));
    parent->connect(Se, OJ2);
    parent->computeStateEquation();
    */
}
