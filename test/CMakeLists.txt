include_directories(BEFORE  
    ${PROJECT_SOURCE_DIR}/src
)

add_executable(
  basic_test
  basic_test.cpp
)
target_link_libraries(
  basic_test
  gtest_main
  ${PROJECT_NAME}
)

include(GoogleTest)
gtest_discover_tests(basic_test)

# executable for testing
add_executable(bgtest test.cpp)
target_link_libraries(bgtest ${PROJECT_NAME})

#add_executable(utest unitstest.cpp)
#target_link_libraries(utest ${PROJECT_NAME})

#add_executable(RLC rlc.cpp)
#target_link_libraries(RLC ${PROJECT_NAME})
#add_test(RLC ${PROJECT_BINARY_DIR}/RLC)

#add_executable(SNF_RLC smithNormalForm.cpp)
#target_link_libraries(SNF_RLC ${PROJECT_NAME})
#add_test(SNF_RLC ${PROJECT_BINARY_DIR}/SNF_RLC)