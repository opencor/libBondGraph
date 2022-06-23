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

#include "bondgraph.hpp"
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace BG;

RCPLIB::RCP<BondGraphInterface> rlc()
{
    newWorkSpace();
    auto ioBondGraph = createBondGraph();

    //Create the resistor
    auto lR = createResistor();
    lR->setParameter("r", "1", "Ohm");
    ioBondGraph->addComponent(lR);

    auto lC1 = createCapacitor();
    lC1->setParameter("C", "1", "Farad");
    lC1->setParameter("q_", "1", "coulomb");
    ioBondGraph->addComponent(lC1);

    //Create the junctions
    //auto lJ0_1 = createZeroJunction();
    auto lJ0_1 = createOneJunction();
    ioBondGraph->addComponent(lJ0_1);

    //Create the bonds
    ioBondGraph->connect(lR, lJ0_1);
    ioBondGraph->connect(lC1, lJ0_1);
    return ioBondGraph;
}

RCPLIB::RCP<BondGraphInterface> multiLink(){
    newWorkSpace();
    auto ioBondGraph = createBondGraph();

    //Create the resistor
    auto lC = createCapacitor();
    ioBondGraph->addComponent(lC);
    nlohmann::json vmap;
    vmap["type"] = "file";
    vmap["compartment"] = "tedt_main";
    vmap["filename"] = "D:/GithubRepositories/bguistandalone/build/Debug/bg/tedt.cellml";
    vmap["target"] = "q_0";          
    vmap["link"] = true;
    //If link is true then importName shoud be provided
    vmap["importName"] = "C_3_CML";
    nlohmann::json tv;
    tv["t"]="time";

    vmap["mapvariables"] = nlohmann::json::array({tv});
    //vmap["mapvariables"] = nlohmann::json::array({"t"});
     //{"compartment":"tedt_main","filename":"tedt.cellml","link":false,"mapvariables":"t","target":"f_3","type":"file"}
    lC->setParameter("i", vmap.dump(), "C");
    ioBondGraph->addComponent(lC);
    return ioBondGraph;
}


RCPLIB::RCP<BondGraphInterface> userdefined()
{
    newWorkSpace();
    auto ioBondGraph = createBondGraph();

    //Create the resistor
    auto lC = createCapacitor();
    
    nlohmann::json vmap;
    vmap["type"] = "file";
    vmap["compartment"] = "A1";
    vmap["filename"] = "C:/Users/rjag008/Desktop/BGtest/Cellmlimportdemo/multiconnectiontest.cellml";
    vmap["target"] = "C";          
    vmap["link"] = false;
    vmap["importName"] = "C_1_CML";
    //If link is true then importName shoud be provided
    nlohmann::json tv;
    tv["t"]="time";

    //Required for usedefined where t is involved
    //vmap["mapvariables"] = nlohmann::json::array({tv});

    lC->setParameter("q_", "1.0", "C");
    ioBondGraph->addComponent(lC);

    auto lFS = createPotentialSource();
     //{"compartment":"tedt_main","filename":"tedt.cellml","link":false,"mapvariables":"t","target":"f_3","type":"file"}
    lFS->setParameter("i", vmap.dump(), "C");
    ioBondGraph->addComponent(lFS);

    auto lTF = createTransformer();
    lTF->setParameter("n","2","");
    ioBondGraph->addComponent(lTF);

    //Create the junctions
    auto lJ0_1 = createZeroJunction();
    //auto lJ0_1 = createOneJunction();
    ioBondGraph->addComponent(lJ0_1);

    //Create the bonds
    ioBondGraph->connect(lC, 0, lTF,0);
    ioBondGraph->connect(lTF,1,lJ0_1);
    ioBondGraph->connect(lJ0_1,lFS);
    return ioBondGraph;
}


RCPLIB::RCP<BondGraphInterface> cellmlParameter()
{
    newWorkSpace();
    auto ioBondGraph = createBondGraph();

    //Create the resistor
    auto lR = createResistor();
    lR->setParameter("r", "1", "Ohm");
    ioBondGraph->addComponent(lR);

    auto lC1 = createCapacitor();
    lC1->setParameter("C", "1", "Farad");
    nlohmann::json vmap;
    vmap["type"] = "file";
    vmap["compartment"] = "A1";
    vmap["filename"] = "C:/Users/rjag008/Desktop/BGtest/Cellmlimportdemo/multiconnectiontest.cellml";
    vmap["target"] = "C";          
    vmap["link"] = false;
    vmap["stateValue"] = true;
    vmap["importName"] = "C_1_CML";

    lC1->setParameter("q_", vmap.dump(), "coulomb");
    ioBondGraph->addComponent(lC1);

    //Create the junctions
    //auto lJ0_1 = createZeroJunction();
    auto lJ0_1 = createOneJunction();
    ioBondGraph->addComponent(lJ0_1);

    //Create the bonds
    ioBondGraph->connect(lR, lJ0_1);
    ioBondGraph->connect(lC1, lJ0_1);

    return ioBondGraph;
}



// RCPLIB::RCP<BondGraphInterface> rlc()
// {
//     auto ioBondGraph = createBondGraph();

//     //Create the storage
//     auto lI = createInductor();
//     ioBondGraph->addComponent(lI);

//     auto lC1 = createCapacitor();
//     ioBondGraph->addComponent(lC1);

//     //Create the resistor
//     auto lR = createResistor();
//     ioBondGraph->addComponent(lR);

//     auto lE = createConstantVoltageSource();
//     ioBondGraph->addComponent(lE);

//     //Create the junctions
//     auto lJ0_1 = createZeroJunction();
//     ioBondGraph->addComponent(lJ0_1);

//     //Create the bonds
//     ioBondGraph->connect(lR, lJ0_1);
//     ioBondGraph->connect(lI, lJ0_1);
//     ioBondGraph->connect(lC1, lJ0_1);
//     ioBondGraph->connect(lE, lJ0_1);

//     /*
//     SymEngine::vec_basic equations;
//     SymEngine::vec_basic bondEquations;
//     SymEngine::vec_basic constraints;

//     std::tie(equations, bondEquations, constraints) = ioBondGraph->computeStateEquation();
//     std::cout << equations << std::endl;
//     std::cout << constraints << std::endl;
//     */
//     ioBondGraph->computeStateEquation();
//     //std::cout<<ioBondGraph->serializeAsCellML("RLC")<<std::endl;
//     return ioBondGraph;
// }



RCPLIB::RCP<BondGraphInterface> reaction()
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
        auto nm = defaultRegistry().getNameMap();
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
        auto nm = defaultRegistry().getNameMap();
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

    RCPLIB::RCP<BondGraphInterface> bg2 = createBondGraph();
    j.get_to(bg2);

    nlohmann::json j2 = bg2;
    std::cout<<" Compare "<< (j==j2) <<std::endl;
    std::cout<<"First\n\n"<<j<<"\nSecond\n\n"<<j2<<std::endl;
    //loadJson(j);
}

#define SETNAME(x)  x->setName(#x);

RCPLIB::RCP<BondGraphInterface> bgtReaction(){
    auto ioBondGraph = createBondGraph();
    auto C_A = createConcentration();
    SETNAME(C_A)
    C_A->setParameter("a_", "1.0", "mol");
    ioBondGraph->addComponent(C_A);

    auto C_B = createConcentration();
    SETNAME(C_B)
    C_B->setParameter("a_", "2.0", "mol");
    ioBondGraph->addComponent(C_B);


    auto C_C = createConcentration();
    SETNAME(C_C)
    C_C->setParameter("a_", "3.0", "mol");
    ioBondGraph->addComponent(C_C);

    auto R1 = createReaction();
    SETNAME(R1)
    ioBondGraph->addComponent(R1);

    auto flow = createOneJunction();
    SETNAME(flow)
    ioBondGraph->addComponent(flow);

    ioBondGraph->connect(C_A,flow);
    ioBondGraph->connect(C_B,flow);
    //ioBondGraph->connect(flow,R1);
    ioBondGraph->connectInverting(flow,R1,0);
    ioBondGraph->connect(R1,1,C_C);

    return ioBondGraph;
}


RCPLIB::RCP<BondGraphInterface> bgtMMReaction(){
    auto ioBondGraph = createBondGraph();
    auto C_E = createConcentration();
    SETNAME(C_E)
    C_E->setParameter("a_", "1.0", "mol");
    ioBondGraph->addComponent(C_E);

    auto C_C = createConcentration();
    SETNAME(C_C)
    C_C->setParameter("a_", "2.0", "mol");
    ioBondGraph->addComponent(C_C);

    auto S_S = createPotentialSource();
    SETNAME(S_S)
    ioBondGraph->addComponent(S_S);
    
    auto S_P = createPotentialSource();
    SETNAME(S_P)
    ioBondGraph->addComponent(S_P);

    auto R1 = createReaction();
    SETNAME(R1)
    ioBondGraph->addComponent(R1);

    auto R2 = createReaction();
    SETNAME(R2)
    ioBondGraph->addComponent(R2);


    auto ej_E = createZeroJunction();
    SETNAME(ej_E)
    ioBondGraph->addComponent(ej_E);

    auto ej_C = createZeroJunction();
    SETNAME(ej_C)
    ioBondGraph->addComponent(ej_C);

    auto fj_1 = createOneJunction();
    SETNAME(fj_1)
    ioBondGraph->addComponent(fj_1);

    auto fj_2 = createOneJunction();
    SETNAME(fj_2)
    ioBondGraph->addComponent(fj_2);

    ioBondGraph->connect(S_S,fj_1) ; // Bond 1
    ioBondGraph->connectInverting(fj_1,R1,0) ; // Bond 2
    ioBondGraph->connect(R1,1,ej_C) ; // Bond 3
    ioBondGraph->connect(ej_C,C_C) ; // Bond 4
    ioBondGraph->connectInverting(ej_C,R2,0) ; // Bond 5
    ioBondGraph->connect(R2,1,fj_2) ; // Bond 6
    ioBondGraph->connect(fj_2,S_P) ; // Bond 7
    ioBondGraph->connect(fj_2,ej_E) ; // Bond 8
    ioBondGraph->connect(ej_E,C_E) ; // Bond 9
    ioBondGraph->connect(ej_E,fj_1) ; // Bond 10

    return ioBondGraph;
}


RCPLIB::RCP<BondGraphInterface> RCSConstraints(){
    auto ioBondGraph = createBondGraph();
    auto lC1 = createCapacitor();
    ioBondGraph->addComponent(lC1);

    //Create the resistor
    auto lR = createResistor();
    ioBondGraph->addComponent(lR);

    auto lE = createConstantCurrentSource();
    ioBondGraph->addComponent(lE);

    //Create the junctions
    auto lJ0_1 = createZeroJunction();
    ioBondGraph->addComponent(lJ0_1);

    //Create the bonds
    ioBondGraph->connect(lR, lJ0_1);
    ioBondGraph->connect(lJ0_1,lC1);
    ioBondGraph->connect(lE, lJ0_1);
    return ioBondGraph;
}


RCPLIB::RCP<BondGraphInterface> gpcrReaction(){
    /* //Bondgraph tools constitutive equations

    dq_m_RG + k0*r19*q_m_RG + k1*r14*q_m_RGL - r19/(k5*k6*q_m_R*q_i_GaGDPbg) - r14/(k0*k13*q_m_RG*q_0_L)
    
    dq_m_RGL + k1*r14*q_m_RGL + k1*r15*q_m_RGL - k2*k8*r15*q_m_RLBG*q_i_GDP - r14/(k0*k13*q_m_RG*q_0_L)
    
    dq_m_RLBG - k1*r15*q_m_RGL + k2*k8*r15*q_m_RLBG*q_i_GDP + k3*r16*q_m_RGLP - r16/(k12*k2*q_i_GTP*q_m_RLBG)
    
    dq_m_RGLP - k10*k11*k4*r17*q_i_Gbg*q_i_GaGTP*q_i_RL + k3*r16*q_m_RGLP + k3*r17*q_m_RGLP - r16/(k12*k2*q_i_GTP*q_m_RLBG)
    
    dq_i_RL + k10*k11*k4*r17*q_i_Gbg*q_i_GaGTP*q_i_RL - k13*k5*r18*q_0_L*q_m_R - k3*r17*q_m_RGLP + k4*r18*q_i_RL
    
    dq_m_R + k0*r19*q_m_RG + k13*k5*r18*q_0_L*q_m_R - k4*r18*q_i_RL - r19/(k5*k6*q_m_R*q_i_GaGDPbg)
    
    dq_i_GaGDPbg + k0*r19*q_m_RG + k6*r20*q_i_GaGDPbg - r19/(k5*k6*q_m_R*q_i_GaGDPbg) - r20/(k10*k7*q_i_Gbg*q_i_GaGDP)
    
    dq_i_GaGDP - k11*r21*q_i_GaGTP + k6*r20*q_i_GaGDPbg + k7*k9*r21*q_i_GaGDP*q_i_P - r20/(k10*k7*q_i_Gbg*q_i_GaGDP)
    
    dq_i_GDP - k1*r15*q_m_RGL + k12*r22*q_i_GTP + k2*k8*r15*q_m_RLBG*q_i_GDP - r22/(k8*k9*q_i_GDP*q_i_P)
    
    dq_i_P - k11*r21*q_i_GaGTP + k12*r22*q_i_GTP + k7*k9*r21*q_i_GaGDP*q_i_P - r22/(k8*k9*q_i_GDP*q_i_P)
    
    dq_i_Gbg - i_Gbg + k10*k11*k4*r17*q_i_Gbg*q_i_GaGTP*q_i_RL - k3*r17*q_m_RGLP + k6*r20*q_i_GaGDPbg - r20/(k10*k7*q_i_Gbg*q_i_GaGDP)
    
    dq_i_GaGTP - i_GaGTP + k10*k11*k4*r17*q_i_Gbg*q_i_GaGTP*q_i_RL + k11*r21*q_i_GaGTP - k3*r17*q_m_RGLP - k7*k9*r21*q_i_GaGDP*q_i_P
    
    dq_i_GTP + k12*r22*q_i_GTP + k3*r16*q_m_RGLP - r22/(k8*k9*q_i_GDP*q_i_P) - r16/(k12*k2*q_i_GTP*q_m_RLBG)
    
    dq_0_L - i_L + k1*r14*q_m_RGL + k13*k5*r18*q_0_L*q_m_R - k4*r18*q_i_RL - r14/(k0*k13*q_m_RG*q_0_L)

    */

    auto ioBondGraph = createBondGraph();

    //Create the storage
    auto q_m_RG = createConcentration();
    SETNAME(q_m_RG)
    q_m_RG->setParameter("a_", "4.18367004218109", "fmol");
    ioBondGraph->addComponent(q_m_RG);

    auto q_m_RGL = createConcentration();
    SETNAME(q_m_RGL)
    q_m_RGL->setParameter("a_", "21.3028137054491", "fmol");
    ioBondGraph->addComponent(q_m_RGL);

    auto q_m_RLBG = createConcentration();
    SETNAME(q_m_RLBG)
    q_m_RLBG->setParameter("a_", "3.47404955471616", "fmol");
    ioBondGraph->addComponent(q_m_RLBG);

    auto q_m_RGLP = createConcentration();
    SETNAME(q_m_RGLP)
    q_m_RGLP->setParameter("a_", "25.5067002874004", "fmol");
    ioBondGraph->addComponent(q_m_RGLP);

    auto q_i_RL = createConcentration();
    SETNAME(q_i_RL)
    q_i_RL->setParameter("a_", "4.62454245115", "fmol");
    ioBondGraph->addComponent(q_i_RL);

    auto q_m_R = createConcentration();
    SETNAME(q_m_R)
    q_m_R->setParameter("a_", "0.908223959103197", "fmol");
    ioBondGraph->addComponent(q_m_R);

    auto q_i_GaGDPbg = createConcentration();
    SETNAME(q_i_GaGDPbg)
    q_i_GaGDPbg->setParameter("a_", "4.60639295980652", "fmol");
    ioBondGraph->addComponent(q_i_GaGDPbg);

    auto q_i_GaGDP = createConcentration();
    SETNAME(q_i_GaGDP)
    q_i_GaGDP->setParameter("a_", "4.97250314034647", "fmol");
    ioBondGraph->addComponent(q_i_GaGDP);

    auto q_i_GDP = createConcentration();
    SETNAME(q_i_GDP)
    q_i_GDP->setParameter("a_", "6.13197134518977", "fmol");
    ioBondGraph->addComponent(q_i_GDP);

    auto q_i_P = createConcentration();
    SETNAME(q_i_P)
    q_i_P->setParameter("a_", "1.19735119297301", "fmol");
    ioBondGraph->addComponent(q_i_P);
   

    auto q_i_Gbg = createConcentration();
    SETNAME(q_i_Gbg)
    q_i_Gbg->setParameter("a_", "0.926373450446675", "fmol");
    ioBondGraph->addComponent(q_i_Gbg);

    auto q_i_GaGTP = createConcentration();
    SETNAME(q_i_GaGTP)
    q_i_GaGTP->setParameter("a_", "7.34207820952638", "fmol");
    ioBondGraph->addComponent(q_i_GaGTP);

    auto q_i_GTP = createConcentration();
    SETNAME(q_i_GTP)
    q_i_GTP->setParameter("a_", "6.13197134518977", "fmol");
    ioBondGraph->addComponent(q_i_GTP);        

    auto q_0_L = createConcentration();
    SETNAME(q_0_L)
    q_0_L->setParameter("a_", "5.09188920026956", "fmol");
    ioBondGraph->addComponent(q_0_L);  

    auto f_0_u_i_Gbg = createConstantFlowSource();
    SETNAME(f_0_u_i_Gbg)
    f_0_u_i_Gbg->setParameter("i", "0", "fmol/s");
    ioBondGraph->addComponent(f_0_u_i_Gbg);  

    auto f_0_u_i_GaGTP = createConstantFlowSource();
    SETNAME(f_0_u_i_GaGTP)
    f_0_u_i_GaGTP->setParameter("i", "0", "fmol/s");
    ioBondGraph->addComponent(f_0_u_i_GaGTP);  

    auto f_0_u_0_L = createConstantFlowSource();
    SETNAME(f_0_u_0_L)
    f_0_u_0_L->setParameter("i", "0", "fmol/s");
    ioBondGraph->addComponent(f_0_u_0_L);  

    //Create the reaction
    auto Re_m_1 = createReaction();
    SETNAME(Re_m_1)
    ioBondGraph->addComponent(Re_m_1);

    auto Re_m_2 = createReaction();
    SETNAME(Re_m_2)
    ioBondGraph->addComponent(Re_m_2);

    auto Re_m_3 = createReaction();
    SETNAME(Re_m_3)
    ioBondGraph->addComponent(Re_m_3);

    auto Re_m_4 = createReaction();
    SETNAME(Re_m_4)
    ioBondGraph->addComponent(Re_m_4);

    auto Re_m_5 = createReaction();
    SETNAME(Re_m_5)
    ioBondGraph->addComponent(Re_m_5);

    auto Re_m_6 = createReaction();
    SETNAME(Re_m_6)
    ioBondGraph->addComponent(Re_m_6);

    auto Re_m_7 = createReaction();
    SETNAME(Re_m_7)
    ioBondGraph->addComponent(Re_m_7);

    auto Re_m_8 = createReaction();
    SETNAME(Re_m_8)
    ioBondGraph->addComponent(Re_m_8);

    auto Re_m_9 = createReaction();
    SETNAME(Re_m_9)
    ioBondGraph->addComponent(Re_m_9);    

    //Create the junctions
    //Onejunction
    auto v_1_m = createOneJunction();
    SETNAME(v_1_m)
    ioBondGraph->addComponent(v_1_m);

    auto v_2_m = createOneJunction();
    SETNAME(v_2_m)
    ioBondGraph->addComponent(v_2_m);

    auto v_3_m = createOneJunction();
    SETNAME(v_3_m)
    ioBondGraph->addComponent(v_3_m);

    auto v_4_m = createOneJunction();
    SETNAME(v_4_m)
    ioBondGraph->addComponent(v_4_m);

    auto v_5_m = createOneJunction();
    SETNAME(v_5_m)
    ioBondGraph->addComponent(v_5_m);

    auto v_6_m = createOneJunction();
    SETNAME(v_6_m)
    ioBondGraph->addComponent(v_6_m);     

    auto v_7_m = createOneJunction();
    SETNAME(v_7_m)
    ioBondGraph->addComponent(v_7_m);

    auto v_8_m = createOneJunction();
    SETNAME(v_8_m)
    ioBondGraph->addComponent(v_8_m);

    auto v_9_m = createOneJunction();
    SETNAME(v_9_m)
    ioBondGraph->addComponent(v_9_m);      

    //Zero junction
    auto u_L_o = createZeroJunction();
    SETNAME(u_L_o)
    ioBondGraph->addComponent(u_L_o);  

    auto u_RG_m = createZeroJunction();
    SETNAME(u_RG_m)
    ioBondGraph->addComponent(u_RG_m);  

    auto u_RGL_m = createZeroJunction();
    SETNAME(u_RGL_m)
    ioBondGraph->addComponent(u_RGL_m); 

    auto u_RLBG_m = createZeroJunction();
    SETNAME(u_RLBG_m)
    ioBondGraph->addComponent(u_RLBG_m); 

    auto u_RGLP_m = createZeroJunction();
    SETNAME(u_RGLP_m)
    ioBondGraph->addComponent(u_RGLP_m);         

    auto u_RL_m = createZeroJunction();
    SETNAME(u_RL_m)
    ioBondGraph->addComponent(u_RL_m);  

    auto u_R_m = createZeroJunction();
    SETNAME(u_R_m)
    ioBondGraph->addComponent(u_R_m);  

    auto u_GaGDPbg_i = createZeroJunction();
    SETNAME(u_GaGDPbg_i)
    ioBondGraph->addComponent(u_GaGDPbg_i);  

    auto u_GaGDP_i = createZeroJunction();
    SETNAME(u_GaGDP_i)
    ioBondGraph->addComponent(u_GaGDP_i);  

    auto u_GDP_i = createZeroJunction();
    SETNAME(u_GDP_i)
    ioBondGraph->addComponent(u_GDP_i);  

    auto u_P_i = createZeroJunction();
    SETNAME(u_P_i)
    ioBondGraph->addComponent(u_P_i);  

    auto u_Gbg_i = createZeroJunction();
    SETNAME(u_Gbg_i)
    ioBondGraph->addComponent(u_Gbg_i);  

    auto u_GaGTP_i = createZeroJunction();
    SETNAME(u_GaGTP_i)
    ioBondGraph->addComponent(u_GaGTP_i);      

    auto u_GTP_i = createZeroJunction();
    SETNAME(u_GTP_i)
    ioBondGraph->addComponent(u_GTP_i);      

    //Create bonds
    //R1
    ioBondGraph->connectInverting(v_1_m, Re_m_1, 0);
    ioBondGraph->connectInverting(Re_m_1, 1, u_RGL_m);

    //R2
    ioBondGraph->connectInverting(u_RGL_m,Re_m_2, 0);
    ioBondGraph->connectInverting(Re_m_2, 1, v_2_m);

    //R3
    ioBondGraph->connectInverting(v_3_m, Re_m_3, 0);
    ioBondGraph->connectInverting(Re_m_3, 1, u_RGLP_m);

    //R4
    ioBondGraph->connectInverting(Re_m_4, 0, u_RGLP_m);
    ioBondGraph->connectInverting(Re_m_4, 1, v_4_m);

    //R5
    ioBondGraph->connectInverting(Re_m_5, 0, u_RL_m);
    ioBondGraph->connectInverting(Re_m_5, 1, v_5_m);   

    //Reaction R6
    ioBondGraph->connectInverting(v_6_m,Re_m_6, 0);
    ioBondGraph->connectInverting(Re_m_6, 1, u_RG_m);

    //Reaction R7
    ioBondGraph->connectInverting(v_7_m, Re_m_7, 0);
    ioBondGraph->connectInverting(Re_m_7, 1, u_GaGDPbg_i);         

    //Reaction R8
    ioBondGraph->connectInverting(Re_m_8, 0, u_GaGTP_i);
    ioBondGraph->connectInverting(Re_m_8, 1, v_8_m);    

    //Reaction R9
    ioBondGraph->connectInverting(v_9_m, Re_m_9, 0);
    ioBondGraph->connectInverting(Re_m_9, 1, u_GTP_i);  

    //Junction connections
    ioBondGraph->connect(u_RG_m,v_1_m);
    ioBondGraph->connect(v_2_m,u_GDP_i);
    ioBondGraph->connect(v_2_m,u_RLBG_m);  
    ioBondGraph->connect(u_RLBG_m,v_3_m);
    ioBondGraph->connect(v_4_m,u_GaGTP_i);
    ioBondGraph->connect(v_4_m,u_Gbg_i);
    ioBondGraph->connect(v_4_m,u_RL_m);
    ioBondGraph->connect(u_Gbg_i,v_7_m);
    ioBondGraph->connect(v_5_m,u_R_m);
    ioBondGraph->connect(u_R_m,v_6_m);
    ioBondGraph->connect(v_5_m,u_L_o);
    ioBondGraph->connect(u_L_o,v_1_m);
    ioBondGraph->connect(u_GaGDPbg_i,v_6_m);
    ioBondGraph->connect(u_GaGDP_i,v_7_m);
    ioBondGraph->connect(v_8_m,u_P_i);
    ioBondGraph->connect(v_8_m,u_GaGDP_i);
    ioBondGraph->connect(u_P_i,v_9_m);
    ioBondGraph->connect(u_GDP_i,v_9_m);
    ioBondGraph->connect(u_GTP_i,v_3_m);

    //Connect concentrations
    ioBondGraph->connect(u_RG_m,q_m_RG);
    ioBondGraph->connect(u_RGL_m,q_m_RGL);
    ioBondGraph->connect(u_RLBG_m,q_m_RLBG);
    ioBondGraph->connect(u_RGLP_m,q_m_RGLP);
    ioBondGraph->connect(u_RL_m,q_i_RL);
    ioBondGraph->connect(u_R_m,q_m_R);
    ioBondGraph->connect(u_GaGDPbg_i,q_i_GaGDPbg);
    ioBondGraph->connect(u_L_o,q_0_L);
    ioBondGraph->connect(u_GaGDP_i,q_i_GaGDP);
    ioBondGraph->connect(u_P_i,q_i_P);
    ioBondGraph->connect(u_GDP_i,q_i_GDP);
    ioBondGraph->connect(u_Gbg_i,q_i_Gbg);
    ioBondGraph->connect(u_GaGTP_i,q_i_GaGTP);
    ioBondGraph->connect(u_GTP_i,q_i_GTP);

    //Add flow sources
    ioBondGraph->connect(f_0_u_0_L,u_L_o);
    ioBondGraph->connect(u_Gbg_i,f_0_u_i_Gbg);
    ioBondGraph->connect(u_GaGTP_i,f_0_u_i_GaGTP);

    return ioBondGraph;
}

//#include <units.hpp>

int main(int argc, char *argv[])
{
    //std::cout<<" In here "<<std::endl;
    // auto sd = getSupportedPhysicalDomainsAndFactoryMethods();
    // std::cout<<sd.dump();
    // checkAdditions();
    // testJSON();
    // portHamiltonian();

    //auto bg = bgtReaction();
    //auto bg = bgtMMReaction();
    //auto bg = rlc();
    //auto bg = gpcrReaction();
    //auto bg = RCSConstraints();
    //auto bg = userdefined();
    auto bg = cellmlParameter();
    auto eqs = bg->computeStateEquation();
    //auto files = getCellML("BGTReaction",bg,eqs);
    //auto files = getCellML("Reaction",bg,eqs);
    auto files = getCellML("CellMLImport",bg,eqs);
    for(auto c : files){
        std::ofstream react;
        react.open("D:/Temp/"+c.first);
        react << c.second;
        react.close();
    }
    //std::cout<<files["Reaction.cellml"]<<std::endl;
    


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
