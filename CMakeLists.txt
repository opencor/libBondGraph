cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME BondGraph)
set(PROJECT_URL https://opencor.ws)
set(_PROJECT_VERSION 1.1)

cmake_policy(SET CMP0048 NEW)

project(${PROJECT_NAME} VERSION ${_PROJECT_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

# Support testing to find libraries in shared_libary mode
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

set(BUILD_TESTING ON CACHE BOOL "Setup google test framework")
set(ENABLE_LOGGING ON CACHE BOOL "Enable logging")

set(PARSE_AND_CREATE_FACTORY_METHODS OFF CACHE BOOL "Parse BG Element descriptions in the resource folder to create faactory methods")

# Build teuchos
add_subdirectory(src/thirdparty/teuchos)

# Include Teuchos headers:
include_directories(BEFORE ${teuchos_SOURCE_DIR})

# Include Teuchos config file (generated per build):
include_directories(BEFORE ${teuchos_BINARY_DIR})

# Inform symengine that we are using Teuchos
add_definitions(-DTHIRDPARTY_TEUCHOS)

# SymEngine Requires boost multiprecision 1.77.0 or higher
set(BG_BOOSTMP_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/thirdparty/boostmp_1_77")

# Link with Symengine
add_subdirectory("${PROJECT_SOURCE_DIR}/src/thirdparty/symengine")

#Eigen library for Port Hamiltonian computation
set(EIGEN_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/thirdparty/eigen/")

# Link with LLNL units for normal build, use patched version for webassembly
if(NOT EMSCRIPTEN)
  
  set(UNITS_INSTALL OFF CACHE INTERNAL "")
  set(UNITS_BUILD_OBJECT_LIBRARY ON CACHE INTERNAL "")
  set(UNITS_BUILD_STATIC_LIBRARY OFF CACHE INTERNAL "")
  set(UNITS_BUILD_SHARED_LIBRARY OFF CACHE INTERNAL "")
  set(UNITS_BUILD_WEBSERVER OFF CACHE INTERNAL "")
  add_subdirectory("${PROJECT_SOURCE_DIR}/src/thirdparty/units")

  # hide units variables
  set(UNITS_HEADER_ONLY OFF CACHE INTERNAL "hide this!" FORCE)
  set(UNITS_BUILD_OBJECT_LIBRARY ${UNITS_BUILD_OBJECT_LIBRARY} CACHE INTERNAL "hide this!" FORCE)
  set(UNITS_NAMESPACE ${UNITS_NAMESPACE} CACHE INTERNAL "hide this!" FORCE)


  # Hide splog variables
  set(BUILD_TESTS OFF CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_ALL ${SPDLOG_BUILD_ALL} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_BENCH ${SPDLOG_BUILD_BENCH} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_EXAMPLE_HO ${SPDLOG_BUILD_EXAMPLE_HO} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_SHARED OFF CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_TESTS ${SPDLOG_BUILD_TESTS} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_TESTS_HO ${SPDLOG_BUILD_TESTS_HO} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_WARNINGS ${SPDLOG_BUILD_WARNINGS} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_CLOCK_COARSE ${SPDLOG_CLOCK_COARSE} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_DISABLE_DEFAULT_LOGGER ${SPDLOG_DISABLE_DEFAULT_LOGGER} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_ENABLE_PCH ${SPDLOG_ENABLE_PCH} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_DISABLE_DEFAULT_LOGGER ${SPDLOG_DISABLE_DEFAULT_LOGGER} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_FMT_EXTERNAL ${SPDLOG_FMT_EXTERNAL} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_FMT_EXTERNAL_HO ${SPDLOG_FMT_EXTERNAL_HO} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_INSTALL ${SPDLOG_INSTALL} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_NO_ATOMIC_LEVELS ${SPDLOG_NO_ATOMIC_LEVELS} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_NO_EXCEPTIONS ${SPDLOG_NO_EXCEPTIONS} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_NO_THREAD_ID ${SPDLOG_NO_THREAD_ID} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_NO_TLS ${SPDLOG_NO_TLS} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_PREVENT_CHILD_FD ${SPDLOG_PREVENT_CHILD_FD} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_TIDY ${SPDLOG_TIDY} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_SANITIZE_ADDRESS ${SPDLOG_SANITIZE_ADDRESS} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_TRACE ${SPDLOG_LEVEL_NAME_TRACE} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_OFF ${SPDLOG_LEVEL_NAME_OFF} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_ERROR ${SPDLOG_LEVEL_NAME_ERROR} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_CRITICAL ${SPDLOG_LEVEL_NAME_CRITICAL} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_DEBUG ${SPDLOG_LEVEL_NAME_DEBUG} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_INFO ${SPDLOG_LEVEL_NAME_INFO} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_LEVEL_NAME_WARNING ${SPDLOG_LEVEL_NAME_WARNING} CACHE INTERNAL "hide this!" FORCE)
  set(SPDLOG_BUILD_EXAMPLE ${SPDLOG_BUILD_EXAMPLE} CACHE INTERNAL "hide this!" FORCE)

  if(BUILD_TESTING)
    # Setup google test
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
    )

    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # So that GoogleTest does not get installed with this project
    option(INSTALL_GTEST "Enable installation of googletest." OFF)
    FetchContent_MakeAvailable(googletest)
  endif()

endif(NOT EMSCRIPTEN)

# Load the Units.cellml file into library
include(cmake/vtkEncodeString.cmake)
set(UNITCELLMLVARNAME cellmlunits)
vtk_encode_string(
  INPUT "src/resource/Units.cellml"
  NAME ${UNITCELLMLVARNAME}
  SOURCE_OUTPUT UNITSCELLMLSOURCE
)

configure_file("${PROJECT_SOURCE_DIR}/src/cellmlunitsmap.h.in" "${PROJECT_SOURCE_DIR}/src/cellmlunitsmap.h")

# Export the metadatajson generated by generateFactoryMethods.py
set(ELEMENTMETADATAJSON metadatastring)
vtk_encode_string(
  INPUT "src/resource/metadata.json"
  NAME ${ELEMENTMETADATAJSON}
  SOURCE_OUTPUT ELEMENTMETADATAJSONSOURCE
)

set(SRC
  src/thirdparty/tinyxml2.cpp
  src/componentregistry.cpp
  src/logging.cpp
  src/utilities.cpp
  src/Port.cpp
  src/Serialisation.cpp
  src/Elementsbase.cpp
  src/Elements.cpp
  src/Bond.cpp
  src/bondgraph.cpp
  src/friends.cpp
  ${UNITSCELLMLSOURCE}
  ${ELEMENTMETADATAJSONSOURCE}
)



# Generate a bondgraph library
if(WIN32)
  add_definitions(-DWIN_EXPORT)
endif()

if(PARSE_AND_CREATE_FACTORY_METHODS)
  find_package(Python3 COMPONENTS Interpreter Development)
  add_custom_target(
    factorymethods ALL
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/resource/generateFactoryMethods.py
    BYPRODUCTS friends.h factorymethods.h factortmethods.cpp
    COMMENT "Generating Domain specific factory methods for bg components"
  )

  #add_dependencies(${PROJECT_NAME} factorymethods)
endif()

if(NOT EMSCRIPTEN)
  include(GNUInstallDirs)
  if(MSVC)
    # Required for PHS build as the size of the obj file is large resulting in C1128
    add_compile_options("/bigobj")
  endif()

  add_library(${PROJECT_NAME} SHARED ${SRC})
  set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty/units/units)
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty/)
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src/thirdparty/)

  if(PARSE_AND_CREATE_FACTORY_METHODS)
    add_dependencies(${PROJECT_NAME} factorymethods)
  endif()

  target_include_directories(${PROJECT_NAME}
    PRIVATE
    ${CMAKE_CURRENT_BUILD_DIR}
    ${BG_BOOSTMP_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
  )

  add_executable(
    phstest
    phstest.cpp
  )

  target_link_libraries(phstest
    teuchos
    symengine
    units-object
    ${PROJECT_NAME}
  )

  target_link_libraries(${PROJECT_NAME} PRIVATE teuchos symengine units-object)

  set_property(TARGET ${PROJECT_NAME} PROPERTY VERSION ${_PROJECT_VERSION})
  set_property(TARGET ${PROJECT_NAME} PROPERTY SOVERSION 1)
  set_property(TARGET ${PROJECT_NAME} PROPERTY
    "INTERFACE_${PROJECT_NAME}_MAJOR_VERSION" 1)
  set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY
    COMPATIBLE_INTERFACE_STRING "${PROJECT_NAME}_MAJOR_VERSION"
  )

  if(ENABLE_LOGGING)
    # Link with spdlog
    add_subdirectory("${PROJECT_SOURCE_DIR}/src/thirdparty/spdlog")
    target_link_libraries(${PROJECT_NAME} PRIVATE "$<BUILD_INTERFACE:spdlog::spdlog_header_only>")
    add_definitions(-DENABLE_LOGGING)
  endif()

  if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
  endif()

  install(TARGETS ${PROJECT_NAME}
    EXPORT BondGraphTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  set(BondGraph_HEADERS
    src/bondgraph.hpp
    src/export.h
    src/RCP.h
    src/factorymethods.h
  )

  install(FILES ${BondGraph_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Headers)
  install(FILES src/thirdparty/json.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/thirdparty)

  install(EXPORT BondGraphTargets
    FILE BondGraphTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
  )
  include(CMakePackageConfigHelpers)
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/BondGraphConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
  )

  # generate the version file for the config file
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/BondGraphConfigVersion.cmake"
    VERSION "${_PROJECT_VERSION}"
    COMPATIBILITY AnyNewerVersion
  )

  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/BondGraphConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/BondGraphConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
  )
else(NOT EMSCRIPTEN)
#https://gist.github.com/WesThorburn/00c47b267a0e8c8431e06b14997778e4
#emcmake cmake ..  -DCMAKE_INSTALL_PREFIX=`pwd`/install
#emmake make
  message(STATUS "Building WASM library.")
  message(STATUS "Run emmake make to build the assembly!")

  set(units_source_files 
    #${PROJECT_SOURCE_DIR}/src/unitproxy/unitconv.cpp 
    ${PROJECT_SOURCE_DIR}/src/unitproxy/units.cpp 
    ${PROJECT_SOURCE_DIR}/src/unitproxy/r20_conv.cpp 
    ${PROJECT_SOURCE_DIR}/src/unitproxy/x12_conv.cpp 
    ${PROJECT_SOURCE_DIR}/src/unitproxy/commodities.cpp 
  )

  add_library(${PROJECT_NAME}_WASM_LIB STATIC 
              ${units_source_files}
              ${SRC})

  set_property(TARGET ${PROJECT_NAME}_WASM_LIB PROPERTY POSITION_INDEPENDENT_CODE ON)

  include_directories(BEFORE ${PROJECT_SOURCE_DIR}/src/unitproxy/)
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty/)
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src/thirdparty/)

  if(PARSE_AND_CREATE_FACTORY_METHODS)
    add_dependencies(${PROJECT_NAME}_WASM_LIB factorymethods)
  endif()


  target_include_directories(${PROJECT_NAME}_WASM_LIB
    PRIVATE
    ${CMAKE_CURRENT_BUILD_DIR}
    ${BG_BOOSTMP_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
  )

  target_link_libraries(${PROJECT_NAME}_WASM_LIB PRIVATE teuchos symengine)


  add_executable(${PROJECT_NAME}_WASM 
               src/wasm.cpp
               )
  
  set(CMAKE_EXECUTABLE_SUFFIX ".js")
  #set_target_properties(${PROJECT_NAME}_WASM PROPERTIES LINK_FLAGS "--no-entry -fsanitize=undefined -sWASM -sALLOW_MEMORY_GROWTH -sNO_DISABLE_EXCEPTION_CATCHING -sMODULARIZE --bind -s EXPORT_NAME='createLibBondGraph'")
  set_target_properties(${PROJECT_NAME}_WASM PROPERTIES LINK_FLAGS "--no-entry -fsanitize=undefined -sWASM -sALLOW_MEMORY_GROWTH -sNO_DISABLE_EXCEPTION_CATCHING --bind")
  target_include_directories(${PROJECT_NAME}_WASM
    PRIVATE
    ${CMAKE_CURRENT_BUILD_DIR}
    ${BG_BOOSTMP_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
  )

  target_link_libraries(${PROJECT_NAME}_WASM ${PROJECT_NAME}_WASM_LIB)

endif(NOT EMSCRIPTEN)