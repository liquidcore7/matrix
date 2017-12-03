cmake_minimum_required(VERSION 3.5)

project(matrixlib_test)
enable_testing()

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

add_executable(matrixlib_test tests.cpp)
# indicates the include paths
target_include_directories(matrixlib_test PRIVATE ${Boost_INCLUDE_DIRS})
# indicates the shared library variant
target_compile_definitions(matrixlib_test PRIVATE "BOOST_TEST_DYN_LINK=1")
# indicates the link paths
target_link_libraries(matrixlib_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# declares a test with our executable
add_test(NAME runTest COMMAND matrixlib_test)