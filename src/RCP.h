#ifndef BONDGRAPH_RCP_H_
#define BONDGRAPH_RCP_H_

#ifdef THIRDPARTY_TEUCHOS
#    include "thirdparty/teuchos/Teuchos_RCP.hpp"
#    define RCPLIB Teuchos
#else
#    include <symengine/symengine_rcp.h>
#    define RCPLIB SymEngine
#endif

#endif