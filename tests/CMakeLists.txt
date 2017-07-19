cmake_minimum_required(VERSION 3.0)
project(unit_tests)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/Modules/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(
        CMAKE_INCLUDE_PATH

        ${CMAKE_SOURCE_DIR}/../dist/include/
        ${CMAKE_SOURCE_DIR}/../dist/include/rdkit/
        ${CMAKE_SOURCE_DIR}/../dist/include/molpher/
        ${CMAKE_INCLUDE_PATH}
)

set(
        CMAKE_LIBRARY_PATH

        ${CMAKE_SOURCE_DIR}/../dist/lib/
        ${CMAKE_LIBRARY_PATH}
)

find_package(CppUnit REQUIRED)
include_directories(${CPPUNIT_INCLUDE_DIR})
set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})

find_package(MolpherLib)
include_directories(${MOLPHER_INCLUDE_DIR})
set(TEST_LIBS ${TEST_LIBS} ${MOLPHER_LIBRARY})

find_package(RDKit REQUIRED)
if(NOT RDKIT_FOUND)
    message(FATAL_ERROR "RDKit library not found")
endif()
include_directories(${RDKIT_INCLUDE_DIR})
set(TEST_LIBS ${TEST_LIBS} ${RDKIT_LIBRARIES})

set(Boost_USE_MULTITHREADED ON)
find_package(
        Boost REQUIRED

        date_time
        filesystem
        program_options
        regex
        wserialization
        serialization
        signals
        thread
        system
)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost library not found")
endif()
include_directories(${Boost_INCLUDE_DIRS})
set(TEST_LIBS ${TEST_LIBS} ${Boost_LIBRARIES})
add_definitions("-DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB")

set(TEST_SRC minimal_test)
set(TEST_HELPERS helpers)
add_executable(
        minimal_test
        ${TEST_SRC}/MinimalTest.cpp
        ${TEST_SRC}/MinimalTest.hpp
        ${TEST_SRC}/MinimalTestRunner.cpp
        ${TEST_HELPERS}/io/stdout.hpp
        ${TEST_HELPERS}/io/stdout.cpp
)
include_directories(${TEST_HELPERS})
target_link_libraries(
        minimal_test
        ${TEST_LIBS}
)