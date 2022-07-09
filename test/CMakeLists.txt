include_directories(BEFORE
  ${PROJECT_SOURCE_DIR}/src
  ${BG_BOOSTMP_INCLUDE_DIR}
)

add_executable(
  basic_test
  basic_test.cpp
)
target_link_libraries(
  basic_test
  gtest_main
  teuchos
  symengine
  units-object
  ${PROJECT_NAME}
)

include(GoogleTest)
gtest_discover_tests(basic_test)

set(DLLPATH "${CMAKE_CURRENT_BINARY_DIR}/src/thirdparty/symengine/Debug/;${CMAKE_CURRENT_BINARY_DIR}/src/thirdparty/teuchos/Debug/;${CMAKE_CURRENT_BINARY_DIR}/src/thirdparty/units/units/Debug/;")

if(UNIX)
  string(REPLACE ";" ":" DLLPATH "${DLLPATH}")
endif()
