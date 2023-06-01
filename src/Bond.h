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

#pragma once

#include "bondgraph.hpp"
#include "Elements.h"
#include <string>


namespace BG {
    
class BondGraph;

/*! \brief %Bond element of a bond graph.
	*  A bond link two Component via Port. 
	*/
class  Bond: public BondInterface
{
private:
    long int mId;
    std::string mName;
    RCPLIB::RCP<PortInterface> mFromPort; //!< Origin port of the bond.
    RCPLIB::RCP<PortInterface> mToPort; //!< Destination port of the bond.
    
    BondStorageType mStorageType; //!< Storage type to complete the group type.
    ComponentGroup mGroup; //!< Bond group association.

    int mStorageIndex; //!< Index identifying the storage component
    
    Bond(long int inId = -1);

public:
    std::string bid;            //! Unique id for the bond
    virtual ~Bond();
    //! Get the incoming bond port.
     RCPLIB::RCP<PortInterface> getFromPort() const;
    //! Get the outcoming bond port.
     RCPLIB::RCP<PortInterface> getToPort() const;
    //! Get the other end of the bond.
     RCPLIB::RCP<PortInterface> getOtherEnd(const RCPLIB::RCP<PortInterface> inPort) const;
    //! Set the identification number of the bond.
     void setId(long int inId);
    //! Get the identification number of the bond.
     long int getId();
    //! Connect two ports together.
     void connect(const RCPLIB::RCP<PortInterface> &inFromPort, const RCPLIB::RCP<PortInterface> &inToPort);
    //! Return the group type of the bond.
    ComponentGroup getGroup() const;
    BondStorageType getStorageType() const;
    //! Set the index of the associated storage component
    void setStorageIndex(int inIndex);
    //! Get the index of the associated storage component
    int getStorageIndex() const;

    friend  bool operator==(const Bond& lhs, const Bond& rhs);
    friend  RCPLIB::RCP<BondInterface> createBond(long int id);
    friend  std::ostream &operator<<(std::ostream &os, const Bond &p);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<BondInterface> &p);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<PortInterface> &p);
    friend void from_json(const nlohmann::json &j, RCPLIB::RCP<BondGraphInterface> &p);
};

   RCPLIB::RCP<BondInterface> createBond(long int id=1);
} // namespace BG

