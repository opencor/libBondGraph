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

#include "Port.h"

namespace BG {
    Port::Port(bool inPower) 
    {
        weight = 1.0;
        receivesPower = inPower;
        portNo = portCounter++;
        id = ""+portNo;
    }

    Port::~Port() {
        
    };

    void Port::setId(std::string id_){
        id = id_;
    }

    std::string Port::getId(){
        return id;
    }

    bool Port::receivePower(){
        return receivesPower;
    }

    //!  Connect to a component.
    void Port::connect(const RCPLIB::RCP<BGElement>  &inComponent)
    {
        mComponent = inComponent;
    }
    //!  Connect to a bond.
    void Port::connectBond(const RCPLIB::RCP<BondInterface> &inBond)
    {
        mBond = inBond;
    }
    //! Release component
    void Port::release(){
        mComponent.reset();
        mBond.reset();
    }

    //Set dofIndex
    void Port::setDofIndex(long int ix)
    {
        mdofIndex = ix;
    }
    //Get dofIndex
    long int Port::dofIndex()
    {
        return mdofIndex;
    }
    //!  Get the connected component.
    RCPLIB::RCP<BGElement>  Port::getComponent() const
    {
        return mComponent;
    }
    //!  Get the connected bond.
    RCPLIB::RCP<BondInterface> Port::getBond() const
    {
        return mBond;
    }
    //Set weight
    void Port::setWeight(double wt)
    {
        weight = wt;
    }
    //! Get weight
    double Port::getWeight() const
    {
        return weight;
    };

     RCPLIB::RCP<PortInterface> createPort(bool inPower){
        return RCPLIB::rcp(new Port(inPower));
    }

    int Port::portCounter = 0;

} // namespace BG
