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

namespace BG {
/*! \brief Bond Graph Port
    *	Store connected components, dof index and weight
    */

class  Port : public PortInterface
{
private:
    RCPLIB::RCP<BGElement>  mComponent; //!< Attached component
    RCPLIB::RCP<BondInterface> mBond; //!< Attached bond
    long int mdofIndex; //!< Index into dof dimensions
    double weight; //!< Port weight -1.0 for inverting
    bool receivesPower; //! True if power arrow is incident on this port
    Port(bool inPower = false);
    
public:
    int portNo = -1;
    std::string id = "";
    ~Port();
    //! Set identification string
     void setId(std::string id);
    //! Get identification string
     std::string getId();

    //! Get Power direction
     bool receivePower();
    //!  Connect to a component.
     void connect(const RCPLIB::RCP<BGElement>  &inComponent);
    //!  Connect to a bond.
     void connectBond(const RCPLIB::RCP<BondInterface> &inBond);
    //! Release
     void release();
    //Set dofIndex
     void setDofIndex(long int ix);
    //Get dofIndex
     long int dofIndex();
    //!  Get the connected component.
     RCPLIB::RCP<BGElement>  getComponent() const;
    //!  Get the connected bond.
     RCPLIB::RCP<BondInterface> getBond() const;
    //Set weight
     void setWeight(double wt);
    //! Get weight
    double getWeight() const;
    friend  RCPLIB::RCP<PortInterface> createPort(bool inPower);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<PortInterface> &p);
    friend RCPLIB::RCP<BGElement>  loadJson(const nlohmann::json &j, const RCPLIB::RCP<BGElementData> &proxy);
private:
    static int portCounter;
};

 RCPLIB::RCP<PortInterface> createPort(bool inPower);
} // namespace BG
