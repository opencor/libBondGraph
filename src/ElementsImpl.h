#ifndef __ELEMENTS_IMPL_H__
#define __ELEMENTS_IMPL_H__
#include "RCP.h"
#include "thirdparty/json.hpp"
#include <symengine/parser.h>
#include <vector>

namespace BG {

/*! \brief Numeric values with units and magnitude information
        *	Structure to store state, parameter values along with the SI unit information
        *   SI Unit in fundamental units is stored as an 7 dimensional array. The array entries determin the power
        *   The units in order are Seconds,Meters,Kilogram,Kelvin,Ampere,Mole,Candela, when all the entries are 0, the quantity is physically dimensionless
        */
class EXPORTED Value
{
public:
    std::string prefix; //Local name, value against which local parameter search will be performed
    std::string name; //Dof dependent name
    std::string value;
    units::precise_unit units;
    bool universalConstant;
    Value(std::string _name, double v = 0.0)
        : prefix(_name)
        , name(_name)
        , universalConstant(false)
    {
        value = std::to_string(v);
        units = units::unit_from_string("dimless");
    }
    Value(std::string _name, std::string expr = "0.0")
        : prefix(_name)
        , name(_name)
        , value(expr)
        , universalConstant(false)
    {
        units = units::unit_from_string("dimless");
    }
    Value(const Value &v)
    {
        prefix = v.prefix;
        name = v.name;
        value = v.value;
        units = v.units;
        universalConstant = v.universalConstant;
    }
    friend std::ostream &operator<<(std::ostream &os, const Value &val);
};

//Concrete implementation of Bondgraph component, proxy components will hold a pointer to an instance of this
struct ElementImpl
{
    std::string mName = ""; //!< Name of the element
    std::string mId = ""; //!< Element's ID
    std::string hamiltonian = ""; //!< Hamltonian energy description (valid if the element is a Port Hamiltonian)
    PassiveType mElementType = notype; //!< Type of element
    ComponentGroup mComponentGroup = nogroup; //!< Component group association.
    PhysicalDomain mDomain = ""; //! Physical Domain
    unsigned int numStates = 0; //! Number of states for this Bondgraph element
    bool allocatedPorts = false; //! True if ports are preallocated
    bool portModified = false; //! True if port assignment was modified, used by junctions
    std::vector<RCPLIB::RCP<Port>> mPorts = {}; //!< Ports connected to the component. Size depends on the type of component.
    std::vector<std::string> constitutiveEq = {}; //!< Constitutive equations
    std::vector<int> constitutiveEqIndex = {}; //!< Index of constitutive equation in connectivity matrix
    std::vector<RCPLIB::RCP<Value>> mParameter = {}; //!Array of state and parameter values associated with this bondgraph element
};

void to_json(nlohmann::json &j, const ElementImpl &p);
void from_json(const nlohmann::json &j, ElementImpl &p);
void to_json(nlohmann::json &j, const RCPLIB::RCP<Port> &p);
void to_json(nlohmann::json &j, const RCPLIB::RCP<Value> &p);

} // namespace BG
#endif