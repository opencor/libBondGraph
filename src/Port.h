#ifndef PORT_H_
#define PORT_H_
#include "Elements.h"
#include "RCP.h"
#include "export.h"
#include "thirdparty/json.hpp"

namespace BG {
/*! \brief Bond Graph Port
    *	Store connected components, dof index and weight
    */

class Bond;
class Reaction;
class BondGraph;

class EXPORTED Port
{
private:
    RCPLIB::RCP<BondGraphElementBase> mComponent; //!< Attached component
    RCPLIB::RCP<Bond> mBond; //!< Attached bond
    long int mdofIndex; //!< Index into dof dimensions
    double weight; //!< Port weight -1.0 for inverting
    bool receivesPower; //! True if power arrow is incident on this port
    Port(bool inPower = false);
    static int portCounter;

public:
    int portNo = -1;
    ~Port();
    //! Get Power direction
    bool receivePower();
    //!  Connect to a component.
    void connect(const RCPLIB::RCP<BondGraphElementBase> &inComponent);
    //!  Connect to a bond.
    void connectBond(const RCPLIB::RCP<Bond> &inBond);
    //! Release
    void release();
    //Set dofIndex
    void setDofIndex(long int ix);
    //Get dofIndex
    long int dofIndex();
    //!  Get the connected component.
    RCPLIB::RCP<BondGraphElementBase> getComponent() const;
    //!  Get the connected bond.
    RCPLIB::RCP<Bond> getBond() const;
    //Set weight
    void setWeight(double wt);
    //! Get weight
    const double getWeight() const;
    friend RCPLIB::RCP<Port> createPort(bool inPower = false);
    friend void to_json(nlohmann::json &j, const RCPLIB::RCP<Port> &p);
    friend RCPLIB::RCP<BondGraphElementBase> loadJson(const nlohmann::json &j, const RCPLIB::RCP<ElementImpl> &proxy);
};

EXPORTED RCPLIB::RCP<Port> createPort(bool inPower);
} // namespace BG

#endif