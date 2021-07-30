#include "Bond.h"
#include "Exceptions.h"
#include <sstream>

namespace BG {

Bond::Bond(long int inId)
{
    bid = std::to_string(rand());
    setId(inId);
}

Bond::~Bond()
{
    mFromPort.reset();
    mToPort.reset();
}

void Bond::setId(long int inId)
{
    mId = inId;
    std::ostringstream lNameStream;
    lNameStream << "BOND" << inId;
    mName = lNameStream.str();
}

void Bond::connect(const RCPLIB::RCP<Port> &inFromPort, const RCPLIB::RCP<Port> &inToPort)
{
    mFromPort = inFromPort;
    mFromPort->connectBond(RCPLIB::rcp(this, false));

    mToPort = inToPort;
    mToPort->connectBond(RCPLIB::rcp(this, false));
}

RCPLIB::RCP<Port> Bond::getOtherEnd(const RCPLIB::RCP<Port> inPort) const
{
    if (mFromPort == inPort)
        return mToPort;
    else if (mToPort == inPort)
        return mFromPort;
    else
        throw BGException("Port not related to this bond");
}

RCPLIB::RCP<Bond> createBond(long int id)
{
    return RCPLIB::rcp(new Bond(id));
}

//! Get the incoming bond port.
RCPLIB::RCP<Port> Bond::getFromPort() const
{
    return mFromPort;
}
//! Get the outcoming bond port.
RCPLIB::RCP<Port> Bond::getToPort() const
{
    return mToPort;
}

//! Get the identification number of the bond.
long int Bond::getId()
{
    return mId;
};

//! Return the group type of the bond.
ComponentGroup Bond::getGroup() const
{
    return mGroup;
}
BondStorageType Bond::getStorageType() const
{
    return mStorageType;
}

//! Set the index of the associated storage component
void Bond::setStorageIndex(int inIndex)
{
    mStorageIndex = inIndex;
}
//! Get the index of the associated storage component
int Bond::getStorageIndex() const
{
    return mStorageIndex;
}

bool operator==(const Bond &lhs, const Bond &rhs)
{
    return lhs.bid == rhs.bid;
}

std::ostream &operator<<(std::ostream &os, const Bond &p)
{
    os << "Bond " << p.mId << std::endl;
    os << " From port " << p.mFromPort->getComponent()->getName() << std::endl;
    os << " To port " << p.mToPort->getComponent()->getName() << std::endl;
    return os;
}

} // namespace BG
