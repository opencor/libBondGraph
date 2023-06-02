# Auckland Bioengineering Institute Bondgraph library

A basic bondgraph library of standard one-port elements, transformers, gyrators and junctions. The library also allows developers to use blackbox elements that wrap a complex, non-physical logic to a Source element with both potential and flow variables. Rudimentary support for Port-Hamiltonian's are also provided. The library serializes the bondgraph as JSON, computable code is serialized as CellML (version 1.1).
Existing bondgraphs can also be added as elements to compose larger systems.

# Usage

```
#include  "bondgraph.h"
#include  "Serialisation.h"
#include  "componentregistry.h"
#include  <sstream>
#include  <string>
#include  <tuple>
#include  <vector>
  
using  namespace  BG;
 

int  main(int  argc, char  *argv[])
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

	auto  eqs = ioBondGraph->computeStateEquation();
	//Cellml code 
	auto  files = getCellML("RLC", *ioBondGraph, eqs);

	return 0;
}
```

## Memory management

The code uses [Teuchos RCP](https://docs.trilinos.org/dev/packages/teuchos/doc/html/classTeuchos_1_1RCP.html) for managing memory 

## Dependencies
The library depends on the following packages/libraries 
[SymEngine](https://github.com/symengine/symengine.git)   for symbolic processing

[Units]( https://github.com/LLNL/units.git)  for checking dimensional consistency and computing units 

[spdlog](https://github.com/gabime/spdlog.git) for logging (can be (dis)enabled at build time)

[nlohmann](https://github.com/nlohmann/json.git) for json related data processing

[NGraph](https://math.nist.gov/~RPozo/ngraph/)  for tracking the connectivity information

## Language
C++ 11

## Developer Notes
### MSVC option
Using the Eigen library with SymEngine to calaculate Port Hamiltonian representation of the bondgraph results in a large library. MSVC requires /bigobj compile option ( else C1128 is thrown. )

### WebAssembly via emscripten
The LLNL/units library causes 'memory out of bounds' issue when built with string handling functions (to_string, unit_from_string). This was due to the very large `base_unit_names` map in units.cpp. This is addressed by removing the definitions associated with the following terms
"::us" | "::nautical" | "::imp" | "::chinese" | "::currency" | "::acre" | "::apothecaries"| "::cup" | "::typographic"| "::laboratory" | "::log" | "::clinical" | "::ft" | "::textile" | "::gal" | "::btu*" | "pound*" | "::troy" | "::oz" | "::av" | "::gold"

A file named `base_units.txt` in the src/unit_proxy directory has been created with the original list of base unit definitions, along with the grep command to generate the shortened list. 

This does not impact biophysical bondgraph definitions, and the use of SI units. The non webassembly library does include all definitions.

A file named `index.html` in src directory is provided for testing, the file along with `BondGraph_WASM.js` and `BondGraph_WASM.wasm` should be made available for loading in a browser.