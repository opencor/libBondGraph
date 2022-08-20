cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME BondGraph)
set(PROJECT_URL https://opencor.ws)
set(_PROJECT_VERSION 1.0)

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

# Inform symengine that we are using Teuchos
add_definitions(-DTHIRDPARTY_TEUCHOS)

# SymEngine Requires boost multiprecision 1.77.0 or higher
set(BG_BOOSTMP_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/thirdparty/boostmp_1_77")

# Link with Symengine
add_subdirectory("${PROJECT_SOURCE_DIR}/src/thirdparty/symengine")

# Link with LLNL units
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
mark_as_advanced(BUILD_TESTS)
set(SPDLOG_BUILD_ALL ${SPDLOG_BUILD_ALL} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_ALL)
set(SPDLOG_BUILD_BENCH ${SPDLOG_BUILD_BENCH} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_BENCH)
set(SPDLOG_BUILD_EXAMPLE_HO ${SPDLOG_BUILD_EXAMPLE_HO} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_EXAMPLE_HO)
set(SPDLOG_BUILD_SHARED OFF CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_SHARED)
set(SPDLOG_BUILD_TESTS ${SPDLOG_BUILD_TESTS} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_TESTS)
set(SPDLOG_BUILD_TESTS_HO ${SPDLOG_BUILD_TESTS_HO} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_TESTS_HO)
set(SPDLOG_BUILD_WARNINGS ${SPDLOG_BUILD_WARNINGS} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_WARNINGS)
set(SPDLOG_CLOCK_COARSE ${SPDLOG_CLOCK_COARSE} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_CLOCK_COARSE)
set(SPDLOG_DISABLE_DEFAULT_LOGGER ${SPDLOG_DISABLE_DEFAULT_LOGGER} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_DISABLE_DEFAULT_LOGGER)
set(SPDLOG_ENABLE_PCH ${SPDLOG_ENABLE_PCH} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_ENABLE_PCH)
set(SPDLOG_DISABLE_DEFAULT_LOGGER ${SPDLOG_DISABLE_DEFAULT_LOGGER} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_DISABLE_DEFAULT_LOGGER)
set(SPDLOG_FMT_EXTERNAL ${SPDLOG_FMT_EXTERNAL} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_FMT_EXTERNAL)
set(SPDLOG_FMT_EXTERNAL_HO ${SPDLOG_FMT_EXTERNAL_HO} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_FMT_EXTERNAL_HO)
set(SPDLOG_INSTALL ${SPDLOG_INSTALL} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_INSTALL)
set(SPDLOG_NO_ATOMIC_LEVELS ${SPDLOG_NO_ATOMIC_LEVELS} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_NO_ATOMIC_LEVELS)
set(SPDLOG_NO_EXCEPTIONS ${SPDLOG_NO_EXCEPTIONS} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_NO_EXCEPTIONS)
set(SPDLOG_NO_THREAD_ID ${SPDLOG_NO_THREAD_ID} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_NO_THREAD_ID)
set(SPDLOG_NO_TLS ${SPDLOG_NO_TLS} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_NO_TLS)
set(SPDLOG_PREVENT_CHILD_FD ${SPDLOG_PREVENT_CHILD_FD} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_PREVENT_CHILD_FD)
set(SPDLOG_TIDY ${SPDLOG_TIDY} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_TIDY)
set(SPDLOG_SANITIZE_ADDRESS ${SPDLOG_SANITIZE_ADDRESS} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_SANITIZE_ADDRESS)
set(SPDLOG_LEVEL_NAME_TRACE ${SPDLOG_LEVEL_NAME_TRACE} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_TRACE)
set(SPDLOG_LEVEL_NAME_OFF ${SPDLOG_LEVEL_NAME_OFF} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_OFF)
set(SPDLOG_LEVEL_NAME_ERROR ${SPDLOG_LEVEL_NAME_ERROR} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_ERROR)
set(SPDLOG_LEVEL_NAME_CRITICAL ${SPDLOG_LEVEL_NAME_CRITICAL} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_CRITICAL)
set(SPDLOG_LEVEL_NAME_DEBUG ${SPDLOG_LEVEL_NAME_DEBUG} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_DEBUG)
set(SPDLOG_LEVEL_NAME_INFO ${SPDLOG_LEVEL_NAME_INFO} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_INFO)
set(SPDLOG_LEVEL_NAME_WARNING ${SPDLOG_LEVEL_NAME_WARNING} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_LEVEL_NAME_WARNING)
set(SPDLOG_BUILD_EXAMPLE ${SPDLOG_BUILD_EXAMPLE} CACHE INTERNAL "hide this!" FORCE)
mark_as_advanced(SPDLOG_BUILD_EXAMPLE)

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

# Load the Units.cellml file into library
include(cmake/vtkEncodeString.cmake)
set(UNITCELLMLVARNAME cellmlunits)
vtk_encode_string(
  INPUT "src/resource/Units.cellml"
  NAME ${UNITCELLMLVARNAME}
  SOURCE_OUTPUT UNITSCELLMLSOURCE
)

configure_file("${PROJECT_SOURCE_DIR}/src/cellmlunitsmap.h.in" "${PROJECT_SOURCE_DIR}/src/cellmlunitsmap.h")

# Configure Teuchos using our CMake options:
configure_file(src/Teuchos_config.h.in Teuchos_config.h)

# Export the metadatajson generated by generateFactoryMethods.py
set(ELEMENTMETADATAJSON metadatastring)
vtk_encode_string(
  INPUT "src/resource/metadata.json"
  NAME ${ELEMENTMETADATAJSON}
  SOURCE_OUTPUT ELEMENTMETADATAJSONSOURCE
)

set(SRC
  src/Teuchos_dyn_cast.cpp
  src/Teuchos_Ptr.cpp
  src/Teuchos_RCPNode.cpp
  src/Teuchos_TestForException.cpp
  src/Teuchos_TypeNameTraits.cpp
  src/Teuchos_stacktrace.cpp
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

include(GNUInstallDirs)

# Generate a bondgraph library
if(WIN32)
  add_definitions(-DWIN_EXPORT)
endif()

add_library(${PROJECT_NAME} SHARED ${SRC})
set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty/units/units)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty/)
include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src/thirdparty/)

target_include_directories(${PROJECT_NAME}
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${BG_BOOSTMP_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC symengine units-object)

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
  target_link_libraries(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:spdlog::spdlog_header_only>")
  add_definitions(-DENABLE_LOGGING)
endif()

if(PARSE_AND_CREATE_FACTORY_METHODS)
  find_package(Python3 COMPONENTS Interpreter Development)
  add_custom_target(
    factorymethods ALL
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/resource/generateFactoryMethods.py
    BYPRODUCTS friends.h factorymethods.h factortmethods.cpp
    COMMENT "Generating Domain specific factory methods for bg components"
  )

  add_dependencies(${PROJECT_NAME} factorymethods)
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
  src/Teuchos_any.hpp
  src/Teuchos_Assert.hpp
  src/Teuchos_ConfigDefs.hpp
  src/Teuchos_ConstTypeTraits.hpp
  src/Teuchos_dyn_cast.hpp
  src/Teuchos_ENull.hpp
  src/Teuchos_Exceptions.hpp
  src/Teuchos_getBaseObjVoidPtr.hpp
  src/Teuchos_map.hpp
  src/Teuchos_NullIteratorTraits.hpp
  src/Teuchos_Ptr.hpp
  src/Teuchos_PtrDecl.hpp
  src/Teuchos_RCP.hpp
  src/Teuchos_RCPDecl.hpp
  src/Teuchos_RCPNode.hpp
  src/Teuchos_stacktrace.hpp
  src/Teuchos_TestForException.hpp
  src/Teuchos_toString.hpp
  src/Teuchos_TypeNameTraits.hpp
  src/Teuchos_RCP.hpp
  ${CMAKE_BINARY_DIR}/Teuchos_config.h
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
