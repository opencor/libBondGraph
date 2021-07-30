#include "Port.h"

namespace BG {
    Port::Port(bool inPower) 
    {
        weight = 1.0;
        receivesPower = inPower;
        portNo = portCounter++;
    }

    Port::~Port() {
        
    };

    bool Port::receivePower(){
        return receivesPower;
    }

    //!  Connect to a component.
    void Port::connect(const RCPLIB::RCP<BondGraphElementBase> &inComponent)
    {
        mComponent = inComponent;
    }
    //!  Connect to a bond.
    void Port::connectBond(const RCPLIB::RCP<Bond> &inBond)
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
    RCPLIB::RCP<BondGraphElementBase> Port::getComponent() const
    {
        return mComponent;
    }
    //!  Get the connected bond.
    RCPLIB::RCP<Bond> Port::getBond() const
    {
        return mBond;
    }
    //Set weight
    void Port::setWeight(double wt)
    {
        weight = wt;
    }
    //! Get weight
    const double Port::getWeight() const
    {
        return weight;
    };

    RCPLIB::RCP<Port> createPort(bool inPower){
        return RCPLIB::rcp(new Port(inPower));
    }

    int Port::portCounter = 0;

} // namespace BG
