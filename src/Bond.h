#ifndef BOND_H_
#define BOND_H_

#include "Elements.h"
#include "Port.h"
#include "RCP.h"
#include "export.h"
#include <string>
#include "thirdparty/json.hpp"

namespace BG {
    
class BondGraph;

/*! \brief %Bond element of a bond graph.
	*  A bond link two Component via Port. 
	*/
class EXPORTED Bond
{
private:
    long int mId;
    std::string mName;
    RCPLIB::RCP<Port> mFromPort; //!< Origin port of the bond.
    RCPLIB::RCP<Port> mToPort; //!< Destination port of the bond.
    
    BondStorageType mStorageType; //!< Storage type to complete the group type.
    ComponentGroup mGroup; //!< Bond group association.

    int mStorageIndex; //!< Index identifying the storage component
    
    Bond(long int inId = -1);

public:
    std::string bid;            //! Unique id for the bond
    virtual ~Bond();
    //! Get the incoming bond port.
    RCPLIB::RCP<Port> getFromPort() const;
    //! Get the outcoming bond port.
    RCPLIB::RCP<Port> getToPort() const;
    //! Get the other end of the bond.
    RCPLIB::RCP<Port> getOtherEnd(const RCPLIB::RCP<Port> inPort) const;
    //! Set the identification number of the bond.
    void setId(long int inId);
    //! Get the identification number of the bond.
    long int getId();
    //! Connect two ports together.
    void connect(const RCPLIB::RCP<Port> &inFromPort, const RCPLIB::RCP<Port> &inToPort);
    //! Return the group type of the bond.
    ComponentGroup getGroup() const;
    BondStorageType getStorageType() const;
    //! Set the index of the associated storage component
    void setStorageIndex(int inIndex);
    //! Get the index of the associated storage component
    int getStorageIndex() const;

    friend bool operator==(const Bond& lhs, const Bond& rhs);
    friend RCPLIB::RCP<Bond> createBond(long int id);
    friend std::ostream &operator<<(std::ostream &os, const Bond &p);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<Bond> &p);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<Port> &p);
    friend void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraph> &p);
};

    EXPORTED RCPLIB::RCP<Bond> createBond(long int id=1);
} // namespace BG

#endif
