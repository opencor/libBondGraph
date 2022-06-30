project(symengine)

set(SYMENGINE_MAJOR_VERSION 0)
set(SYMENGINE_MINOR_VERSION 9)
set(SYMENGINE_PATCH_VERSION 0)
set(SYMENGINE_VERSION ${SYMENGINE_MAJOR_VERSION}.${SYMENGINE_MINOR_VERSION}.${SYMENGINE_PATCH_VERSION})

set(INTEGER_CLASS "boostmp" CACHE STRING "Integer class for ${PROJECT_NAME}. ")
set(SYMENGINE_INTEGER_CLASS "BOOSTMP")
set(WITH_BOOST yes)
set(HAVE_SYMENGINE_GMP no)

add_subdirectory(utilities/matchpycpp)

set(SRC
    add.cpp
    basic.cpp
    complex.cpp
    complex_double.cpp
    constants.cpp
    cse.cpp
    cwrapper.cpp
    dense_matrix.cpp
    derivative.cpp
    dict.cpp
    diophantine.cpp
    eval.cpp
    eval_double.cpp
    expand.cpp
    expression.cpp
    fields.cpp
    finitediff.cpp
    functions.cpp
    infinity.cpp
    integer.cpp
    logic.cpp
    matrix.cpp
    monomials.cpp
    mul.cpp
    nan.cpp
    ntheory.cpp
    ntheory_funcs.cpp
    number.cpp
    numer_denom.cpp
    parser/parser_old.cpp
    parser/parser.tab.cc
    parser/parser.cpp
    parser/tokenizer.cpp
    parser/sbml/sbml_parser.tab.cc
    parser/sbml/sbml_parser.cpp
    parser/sbml/sbml_tokenizer.cpp
    polys/basic_conversions.cpp
    polys/msymenginepoly.cpp
    polys/uexprpoly.cpp
    polys/uintpoly.cpp
    polys/uratpoly.cpp
    pow.cpp
    prime_sieve.cpp
    printers/codegen.cpp
    printers/mathml.cpp
    printers/sbml.cpp
    printers/strprinter.cpp
    printers/latex.cpp
    printers/unicode.cpp
    printers/stringbox.cpp
    rational.cpp
    real_double.cpp
    rewrite.cpp
    rings.cpp
    series.cpp
    series_generic.cpp
    sets.cpp
    set_funcs.cpp
    solve.cpp
    sparse_matrix.cpp
    symbol.cpp
    symengine_rcp.cpp
    visitor.cpp
    test_visitors.cpp
    assumptions.cpp
    refine.cpp
    simplify.cpp
    tuple.cpp
)

# Using boost multirecision
set(SRC ${SRC} mp_boost.cpp)

if(WITH_LLVM)
    set(SRC ${SRC} llvm_double.cpp)
endif()

# Needed for "make install"
set(HEADERS
    add.h
    basic.h
    basic-inl.h
    basic-methods.inc
    complex_double.h
    complex.h
    complex_mpc.h
    constants.h
    cwrapper.h
    derivative.h
    dict.h
    diophantine.h
    eval_arb.h
    eval_double.h
    eval.h
    eval_mpc.h
    eval_mpfr.h
    expression.h
    fields.h
    finitediff.h
    flint_wrapper.h
    functions.h
    infinity.h
    integer.h
    lambda_double.h
    llvm_double.h
    logic.h
    matrix.h
    monomials.h
    mp_class.h
    mp_wrapper.h
    mul.h
    nan.h
    ntheory.h
    ntheory_funcs.h
    number.h
    parser.h
    parser/parser.h
    parser/tokenizer.h
    parser/sbml/sbml_parser.h
    parser/sbml/sbml_tokenizer.h
    polys/basic_conversions.h
    polys/cancel.h
    polys/uexprpoly.h
    polys/uintpoly_flint.h
    polys/uintpoly.h
    polys/uintpoly_piranha.h
    polys/upolybase.h
    polys/uratpoly.h
    polys/usymenginepoly.h
    polys/msymenginepoly.h
    pow.h
    prime_sieve.h
    printers/codegen.h
    printers/mathml.h
    printers/sbml.h
    printers/strprinter.h
    printers/latex.h
    printers/unicode.h
    printers.h
    rational.h
    real_double.h
    real_mpfr.h
    rings.h
    serialize-cereal.h
    series_flint.h
    series_generic.h
    series.h
    series_piranha.h
    series_visitor.h
    sets.h
    solve.h
    subs.h
    symbol.h
    symengine_assert.h
    symengine_casts.h
    symengine_exception.h
    symengine_rcp.h
    type_codes.inc
    visitor.h
    test_visitors.h
    assumptions.h
    refine.h
    simplify.h
    utilities/stream_fmt.h
    tuple.h
)

file(GLOB_RECURSE CEREAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "utilities/cereal/include/*")
set(HEADERS ${HEADERS} ${CEREAL_HEADERS})

# Boost.multiprecision build segfaults if -O3 optimization is on
# which it is by default for cmake Release mode. Remove -O3 from flags.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    string(REPLACE "-O3" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

include(GNUInstallDirs)

# Configure SymEngine using our CMake options:
# set(SYMENGINE_CLING_LIBRARY_DIR "\"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\"")
configure_file(symengine_config.h.in symengine_config.h)
configure_file(symengine_config_cling.h.in symengine_config_cling.h)

# Include the config file:
include_directories(BEFORE ${symengine_BINARY_DIR})

# Include the source directory
include_directories(BEFORE ${symengine_SOURCE_DIR})
include_directories(BEFORE ${symengine_BINARY_DIR}/..)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/..)
include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/../teuchos)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../teuchos)
include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/utilities/cereal/include)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/utilities/cereal/include)

add_library(${PROJECT_NAME} STATIC ${SRC})

target_include_directories(${PROJECT_NAME}
    PRIVATE
    ${Boost_INCLUDE_DIR}
)

include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_DEFINITIONS "symengine_EXPORTS")

# if (BUILD_SHARED_LIBS)
# set_target_properties(
# ${PROJECT_NAME}
# PROPERTIES
# VERSION ${SYMENGINE_VERSION}
# SOVERSION ${SYMENGINE_MAJOR_VERSION}.${SYMENGINE_MINOR_VERSION}
# )
# if (WIN32)
# set_target_properties(
# ${PROJECT_NAME}
# PROPERTIES
# RUNTIME_OUTPUT_NAME ${PROJECT_NAME}-${SYMENGINE_MAJOR_VERSION}.${SYMENGINE_MINOR_VERSION}
# )
# endif()
# endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-Wconversion -Wno-sign-conversion")
    set_source_files_properties(parser/parser.tab.cc PROPERTIES COMPILE_FLAGS "-Wno-conversion")
endif()

if(WITH_LLVM AND NOT ${LLVM_VERSION} VERSION_LESS "10.0.0")
    set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF CXX_STANDARD_REQUIRED yes)
endif()

target_link_libraries(${PROJECT_NAME} teuchos)

set_target_properties(${PROJECT_NAME} PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)

if(WITH_LLVM AND NOT SYMENGINE_LLVM_LINK_DOWNSTREAM)
    if(${EXCLUDE_LIBS})
        set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS "-Wl,--exclude-libs,ALL")
    elseif(${UNEXPORTED_SYMBOL})
        set(LLVM_SYMBOLS_AVOID "${LLVM_SYMBOLS_AVOID} -Wl,-unexported_symbol,__ZTVN4llvm* -Wl,-unexported_symbol,__ZNK4llvm*")
        set(LLVM_SYMBOLS_AVOID "${LLVM_SYMBOLS_AVOID} -Wl,-unexported_symbol,__ZN4llvm* -Wl,-unexported_symbol,__ZTIN4llvm*")
        set(LLVM_SYMBOLS_AVOID "${LLVM_SYMBOLS_AVOID} -Wl,-unexported_symbol,_LLVM* -Wl,-unexported_symbol,__ZTSN4llvm*")
        set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS "${LLVM_SYMBOLS_AVOID}")
    endif()
endif()

if(BUILD_SYMENGINE_TESTING)
    add_subdirectory(utilities/catch)
    add_subdirectory(tests)
endif()

# install(TARGETS symengine
# EXPORT BondGraphTargets
# RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
# ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
# LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
# )

# install(FILES ${HEADERS} "${symengine_BINARY_DIR}/symengine_config.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thirdparty/symengine)
