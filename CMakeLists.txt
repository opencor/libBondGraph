cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME BondGraph)
set(PROJECT_URL https://opencor.ws)
set(_PROJECT_VERSION 1.0)

cmake_policy(SET CMP0048 NEW)

project(${PROJECT_NAME} VERSION ${_PROJECT_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(BUILD_TESTING OFF CACHE BOOL "Setup google test framework")
set(ENABLE_LOGGING ON CACHE BOOL "Enable logging")

set(PARSE_AND_CREATE_FACTORY_METHODS OFF CACHE BOOL "Parse BG Element descriptions in the resource folder to create faactory methods")

#Link with Symengine
find_package(SymEngine CONFIG PATHS /hpc/rjag008/Systems_Modelling/usr/lib/cmake/)

#Link with LLNL units
#find_package(units REQUIRED CONFIG PATHS /hpc/rjag008/Systems_Modelling/usr/lib64/cmake/)
set(UNITS_INSTALL OFF CACHE INTERNAL "")
set(UNITS_BUILD_OBJECT_LIBRARY OFF CACHE INTERNAL "")
set(UNITS_BUILD_STATIC_LIBRARY ON CACHE INTERNAL "")
set(UNITS_BUILD_SHARED_LIBRARY OFF CACHE INTERNAL "")
set(UNITS_BUILD_WEBSERVER OFF CACHE INTERNAL "")

add_subdirectory("${PROJECT_SOURCE_DIR}/src/thirdparty/units")


#hide units variables 
set(UNITS_HEADER_ONLY ${UNITS_HEADER_ONLY} CACHE INTERNAL "hide this!" FORCE )
set(UNITS_BUILD_OBJECT_LIBRARY ${UNITS_BUILD_OBJECT_LIBRARY} CACHE INTERNAL "hide this!" FORCE )
set(UNITS_NAMESPACE ${UNITS_NAMESPACE} CACHE INTERNAL "hide this!" FORCE )

set_target_properties(units-static PROPERTIES FOLDER Extern)

#Hide splog variables
set(BUILD_TESTS OFF CACHE INTERNAL "hide this!" FORCE)
set(SPDLOG_BUILD_ALL ${SPDLOG_BUILD_ALL} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_BENCH ${SPDLOG_BUILD_BENCH} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_EXAMPLE_HO ${SPDLOG_BUILD_EXAMPLE_HO} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_SHARED ${SPDLOG_BUILD_SHARED} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_TESTS ${SPDLOG_BUILD_TESTS} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_TESTS_HO ${SPDLOG_BUILD_TESTS_HO} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_WARNINGS ${SPDLOG_BUILD_WARNINGS} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_CLOCK_COARSE ${SPDLOG_CLOCK_COARSE} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_DISABLE_DEFAULT_LOGGER ${SPDLOG_DISABLE_DEFAULT_LOGGER} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_ENABLE_PCH ${SPDLOG_ENABLE_PCH} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_DISABLE_DEFAULT_LOGGER ${SPDLOG_DISABLE_DEFAULT_LOGGER} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_FMT_EXTERNAL ${SPDLOG_FMT_EXTERNAL} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_FMT_EXTERNAL_HO ${SPDLOG_FMT_EXTERNAL_HO} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_INSTALL ${SPDLOG_INSTALL} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_NO_ATOMIC_LEVELS ${SPDLOG_NO_ATOMIC_LEVELS} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_NO_EXCEPTIONS ${SPDLOG_NO_EXCEPTIONS} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_NO_THREAD_ID ${SPDLOG_NO_THREAD_ID} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_NO_TLS ${SPDLOG_NO_TLS} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_PREVENT_CHILD_FD ${SPDLOG_PREVENT_CHILD_FD} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_TIDY ${SPDLOG_TIDY} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_SANITIZE_ADDRESS ${SPDLOG_SANITIZE_ADDRESS} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_TRACE ${SPDLOG_LEVEL_NAME_TRACE} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_OFF ${SPDLOG_LEVEL_NAME_OFF} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_ERROR ${SPDLOG_LEVEL_NAME_ERROR} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_CRITICAL ${SPDLOG_LEVEL_NAME_CRITICAL} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_DEBUG ${SPDLOG_LEVEL_NAME_DEBUG} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_INFO ${SPDLOG_LEVEL_NAME_INFO} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_LEVEL_NAME_WARNING ${SPDLOG_LEVEL_NAME_WARNING} CACHE INTERNAL "hide this!" FORCE )
set(SPDLOG_BUILD_EXAMPLE ${SPDLOG_BUILD_EXAMPLE} CACHE INTERNAL "hide this!" FORCE )

set(SYMENGINE_SKIP_DEPENDENCIES ${SYMENGINE_SKIP_DEPENDENCIES} CACHE INTERNAL "hide this!" FORCE )

if (BUILD_TESTING)
  #Setup google test
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

#Load the Units.cellml file into library
set(UNITSCELLMLFILE "resource/Units.cellml" CACHE INTERNAL "hide this!" FORCE)
file(READ "${PROJECT_SOURCE_DIR}/src/${UNITSCELLMLFILE}" Units_cellml HEX)
configure_file("${PROJECT_SOURCE_DIR}/src/cellmlunitsmap.h.in" "${PROJECT_SOURCE_DIR}/src/cellmlunitsmap.h")

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
)
#Generate a static bondgraph library
add_library(${PROJECT_NAME} STATIC ${SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${SYMENGINE_INCLUDE_DIRS} ${units_INCLUDE_DIRS})
#https://gitlab.kitware.com/cmake/cmake/-/issues/15415#note_634114
target_link_libraries(${PROJECT_NAME} PUBLIC ${SYMENGINE_LIBRARIES} "$<BUILD_INTERFACE:units::units>")

if (ENABLE_LOGGING)
  #Link with spdlog
  add_subdirectory("${PROJECT_SOURCE_DIR}/src/thirdparty/spdlog")
  target_link_libraries(${PROJECT_NAME} PRIVATE "$<BUILD_INTERFACE:spdlog::spdlog_header_only>")
  add_definitions(-DENABLE_LOGGING)
endif()

if (PARSE_AND_CREATE_FACTORY_METHODS)
  find_package (Python3 COMPONENTS Interpreter Development)
  add_custom_target(
    factorymethods ALL
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/resource/generateFactoryMethods.py
    BYPRODUCTS friends.h factorymethods.h factortmethods.cpp
    COMMENT "Generating Domain specific factory methods for bg components"
  )

  add_dependencies(${PROJECT_NAME} factorymethods)
endif()


if (NOT TARGET teuchos)
  add_subdirectory(src/thirdparty/teuchos)
  # Include Teuchos headers:
  include_directories(BEFORE ${teuchos_SOURCE_DIR})
  # Include Teuchos config file (generated per build):
  include_directories(BEFORE ${teuchos_BINARY_DIR})
  target_link_libraries(${PROJECT_NAME} PUBLIC teuchos)
  #Define the flag
  add_definitions(-DTHIRDPARTY_TEUCHOS )
endif()


install(TARGETS ${PROJECT_NAME} 
        EXPORT BondGraphTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}        
        )
install(FILES src/bondgraph.h DESTINATION  ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT BondGraphTargets 
  FILE BondGraphTargets.cmake
  DESTINATION lib/cmake/BondGraph
)

if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()