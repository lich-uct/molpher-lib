cmake_minimum_required(VERSION 3.0)
project(unit_tests)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

find_package(CppUnit REQUIRED)
include_directories(${CPPUNIT_INCLUDE_DIR})
set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})

set(MOLPHER_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../dist/)
find_package(MolpherLib)
include_directories(${MOLPHER_INCLUDE_DIR})
set(TEST_LIBS ${TEST_LIBS} ${MOLPHER_LIBRARY})

add_executable(minimal_test minimal_test/MinimalTest.cpp minimal_test/MinimalTest.hpp minimal_test/MinimalTestRunner.cpp)
target_link_libraries(
        minimal_test
        ${TEST_LIBS}
)