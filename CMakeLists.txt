cmake_minimum_required(VERSION 3.9)
project(molpher_lib)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" "${CMAKE_SOURCE_DIR}/tests/cmake/Modules/")
message("Build type: ${CMAKE_BUILD_TYPE}")

# versioning
file(STRINGS "BUILD.TXT" BUILD_NUMBER LIMIT_COUNT 1)
file(STRINGS "VERSION.TXT" VERSION_NUMBER LIMIT_COUNT 1)
string(REGEX MATCH "^[0-9]+\\.[0-9]+" VERSION_NUMBER_SO "${VERSION_NUMBER}")

# Tool Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")
add_definitions("-D_GLIBCXX_USE_CXX11_ABI=0")
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/.")

# Dependencies

## search paths
set(CONDA_PREFIX $ENV{CONDA_PREFIX})
if(CONDA_PREFIX)
    set(
            CMAKE_PREFIX_PATH

            $ENV{CONDA_PREFIX}
    )
else()
    set(DEPS_DIR "deps")
    set(
            CMAKE_INCLUDE_PATH

            ${DEPS_DIR}/boost/
            ${DEPS_DIR}/rdkit/Code/
            ${DEPS_DIR}/tbb/include/
    )
    set(
            CMAKE_LIBRARY_PATH

            ${DEPS_DIR}/boost/stage/lib/
            ${DEPS_DIR}/rdkit/lib/
            ${DEPS_DIR}/tbb/lib/intel64/gcc4.4/
    )
endif()
message("CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
message("CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
message("CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH}")

## Threads
find_package(Threads REQUIRED)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "Threads library not found")
endif()
include_directories(${Threads_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Threads_LIBRARIES})
add_definitions(${Threads_DEFINITIONS})

## Boost
#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_USE_MULTITHREADED ON)
find_package(
        Boost REQUIRED

#        date_time
        filesystem
#        program_options
#        regex
#        wserialization
        serialization
#        signals
        thread
#        system
)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost library not found")
endif()
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})
add_definitions(${Boost_DEFINITIONS})
add_definitions("-DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB")
message(STATUS "Boost include directories: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries found: ${Boost_LIBRARIES}")

## RDKit
#set(RDKIT_LINK_STATIC ON)
find_package(RDKit REQUIRED)
if(NOT RDKIT_FOUND)
    message(FATAL_ERROR "RDKit library or include files not found")
endif()
include_directories(${RDKIT_INCLUDE_DIR})
set(LIBS ${LIBS} ${RDKIT_LIBRARIES})
add_definitions(${RDKIT_DEFINITIONS})

## TBB
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TBB_USE_DEBUG_BUILD ON)
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(TBB_USE_DEBUG_BUILD OFF)
endif()
message(TBB_USE_DEBUG_BUILD="${TBB_USE_DEBUG_BUILD}")
find_package(TBB REQUIRED)
if(NOT TBB_FOUND)
    message(FATAL_ERROR "TBB library not found")
endif()
include_directories(${TBB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${TBB_LIBRARIES})
add_definitions(${TBB_DEFINITIONS})
message(STATUS "Found TBB include directory: ${TBB_INCLUDE_DIRS}")
message(STATUS "Found TBB libraries at ${TBB_LIBRARIES}")

## CPPUNIT
if(INCLUDE_TESTS)
    find_package(CppUnit)
endif()
if (CPPUNIT_LIBRARIES)
    set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
    message(STATUS "Found CPPUnit include directory: ${CPPUNIT_INCLUDE_DIR}")
    message(STATUS "Found CPPUnit libraries at ${CPPUNIT_LIBRARIES}")
endif ()

# Include Directories
include_directories(
        include/
        src/
)

# Main Target
set(SOURCE_FILES
        include/data_structs/ExplorationData.hpp include/data_structs/ExplorationTree.hpp include/data_structs/MolpherMol.hpp include/operations/callbacks/TraverseCallback.hpp include/operations/ExtendTreeOper.hpp include/operations/FilterMorphsOper.hpp include/operations/FindLeavesOper.hpp include/operations/GenerateMorphsOper.hpp include/operations/PruneTreeOper.hpp include/operations/SortMorphsOper.hpp include/operations/TraverseOper.hpp include/operations/TreeOperation.hpp include/selectors/chemoper_selectors.h include/selectors/fingerprint_selectors.h include/selectors/simcoeff_selectors.h include/SAScore_data_loader.hpp src/core/API/operations/callbacks/TraverseCallback.cpp src/core/API/operations/callbacks/TraverseCallbackImpl.hpp src/core/API/operations/ExtendTreeOper.cpp src/core/API/operations/ExtendTreeOperImpl.hpp src/core/API/operations/FilterMorphsOper.cpp src/core/API/operations/FilterMorphsOperImpl.hpp src/core/API/operations/FindLeavesOper.cpp src/core/API/operations/FindLeavesOperImpl.hpp src/core/API/operations/GenerateMorphsOper.cpp src/core/API/operations/GenerateMorphsOperImpl.hpp src/core/API/operations/PruneTreeOper.cpp src/core/API/operations/PruneTreeOperImpl.hpp src/core/API/operations/SortMorphsOper.cpp src/core/API/operations/SortMorphsOperImpl.hpp src/core/API/operations/TraverseOper.cpp src/core/API/operations/TraverseOperImpl.hpp src/core/API/operations/TreeOperation.cpp src/core/API/operations/TreeOperationImpl.hpp src/core/API/ExplorationData.cpp src/core/API/ExplorationTree.cpp src/core/API/ExplorationTreeImpl.h src/core/API/MolpherMol.cpp src/core/API/MolpherMolImpl.hpp src/core/chem/fingerprintStrategy/AtomPairsFngpr.cpp src/core/chem/fingerprintStrategy/AtomPairsFngpr.hpp src/core/chem/fingerprintStrategy/FingerprintStrategy.cpp src/core/chem/fingerprintStrategy/FingerprintStrategy.h src/core/chem/fingerprintStrategy/MorganFngpr.cpp src/core/chem/fingerprintStrategy/MorganFngpr.hpp src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.hpp src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.hpp src/core/chem/fingerprintStrategy/TopolSingleFngpr.cpp src/core/chem/fingerprintStrategy/TopolSingleFngpr.hpp src/core/chem/fingerprintStrategy/TopolTorsFngpr.cpp src/core/chem/fingerprintStrategy/TopolTorsFngpr.hpp src/core/chem/fingerprintStrategy/VectorFpFngpr.cpp src/core/chem/fingerprintStrategy/VectorFpFngpr.hpp src/core/chem/morphing/CalculateDistances.hpp src/core/chem/morphing/CalculateMorphs.hpp src/core/chem/morphing/Morphing.cpp src/core/chem/morphing/Morphing.hpp src/core/chem/morphing/MorphingData.cpp src/core/chem/morphing/MorphingData.h src/core/chem/morphing/MorphingFtors.cpp src/core/chem/morphing/MorphingFtors.hpp src/core/chem/morphing/ReturnResults.hpp src/core/chem/morphingStrategy/MorphingStrategy.h src/core/chem/morphingStrategy/OpAddAtom.cpp src/core/chem/morphingStrategy/OpAddAtom.hpp src/core/chem/morphingStrategy/OpAddBond.cpp src/core/chem/morphingStrategy/OpAddBond.hpp src/core/chem/morphingStrategy/OpBondContraction.cpp src/core/chem/morphingStrategy/OpBondContraction.hpp src/core/chem/morphingStrategy/OpBondReroute.cpp src/core/chem/morphingStrategy/OpBondReroute.hpp src/core/chem/morphingStrategy/OpInterlayAtom.cpp src/core/chem/morphingStrategy/OpInterlayAtom.hpp src/core/chem/morphingStrategy/OpMutateAtom.cpp src/core/chem/morphingStrategy/OpMutateAtom.hpp src/core/chem/morphingStrategy/OpRemoveAtom.cpp src/core/chem/morphingStrategy/OpRemoveAtom.hpp src/core/chem/morphingStrategy/OpRemoveBond.cpp src/core/chem/morphingStrategy/OpRemoveBond.hpp src/core/chem/simCoefStrategy/AllBitSimCoef.cpp src/core/chem/simCoefStrategy/AllBitSimCoef.hpp src/core/chem/simCoefStrategy/AsymmetricSimCoef.cpp src/core/chem/simCoefStrategy/AsymmetricSimCoef.hpp src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.hpp src/core/chem/simCoefStrategy/CosineSimCoef.cpp src/core/chem/simCoefStrategy/CosineSimCoef.hpp src/core/chem/simCoefStrategy/DiceSimCoef.cpp src/core/chem/simCoefStrategy/DiceSimCoef.hpp src/core/chem/simCoefStrategy/KulczynskiSimCoef.cpp src/core/chem/simCoefStrategy/KulczynskiSimCoef.hpp src/core/chem/simCoefStrategy/McConnaugheySimCoef.cpp src/core/chem/simCoefStrategy/McConnaugheySimCoef.hpp src/core/chem/simCoefStrategy/OnBitSimCoef.cpp src/core/chem/simCoefStrategy/OnBitSimCoef.hpp src/core/chem/simCoefStrategy/RusselSimCoef.cpp src/core/chem/simCoefStrategy/RusselSimCoef.hpp src/core/chem/simCoefStrategy/SimCoefStrategy.h src/core/chem/simCoefStrategy/SokalSimCoef.cpp src/core/chem/simCoefStrategy/SokalSimCoef.hpp src/core/chem/simCoefStrategy/TanimotoSimCoef.cpp src/core/chem/simCoefStrategy/TanimotoSimCoef.hpp src/core/chem/simCoefStrategy/TverskySimCoef.cpp src/core/chem/simCoefStrategy/TverskySimCoef.hpp src/core/chem/ChemicalAuxiliary.cpp src/core/chem/ChemicalAuxiliary.h src/core/chem/SimCoefCalculator.cpp src/core/chem/SimCoefCalculator.hpp src/core/data_structs/ExplorationDataImpl.hpp src/core/data_structs/MolpherMolData.hpp src/core/data_structs/MolpherParam.h src/core/misc/selectors/chemoper_selectors.cpp src/core/misc/selectors/fingerprint_selectors.cpp src/core/misc/selectors/simcoeff_selectors.cpp src/core/misc/global_types.h src/core/misc/inout.cpp src/core/misc/inout.h src/core/misc/iteration_serializer.cpp src/core/misc/iteration_serializer.hpp src/core/misc/SAScore.cpp src/core/misc/SAScore.h src/core/misc/SAScore_data_loader.cpp src/core/misc/SynchRand.cpp src/core/misc/SynchRand.h include/random_seed.hpp src/core/misc/random_seed.cpp include/operations/callbacks/SortMorphsCallback.hpp src/core/API/operations/callbacks/SortMorphsCallbackImpl.hpp src/core/API/operations/callbacks/SortMorphsCallback.cpp include/operations/CleanMorphsOper.hpp src/core/API/operations/CleanMorphsOperImpl.hpp src/core/API/operations/CleanMorphsOper.cpp include/data_structs/MolpherAtom.hpp src/core/API/MolpherAtomImpl.hpp src/core/API/MolpherAtom.cpp include/morphing/operators/MorphingOperator.hpp src/core/morphing/operators/MorphingOperatorImpl.hpp src/core/morphing/operators/MorphingOperator.cpp include/morphing/operators/AddAtom.hpp src/core/morphing/operators/AddAtomImpl.hpp src/core/morphing/operators/AddAtom.cpp include/morphing/AtomLibrary.hpp src/core/morphing/AtomLibrary.cpp src/core/morphing/AtomLibraryImpl.hpp include/morphing/operators/RemoveAtom.hpp src/core/morphing/operators/RemoveAtom.cpp src/core/morphing/operators/RemoveAtomImpl.hpp src/core/misc/utils.hpp src/core/misc/utils.cpp include/morphing/Molpher.hpp src/core/morphing/Molpher.cpp src/core/morphing/MolpherImpl.hpp src/core/morphing/MorphCalculator.cpp src/core/morphing/MorphCalculator.hpp include/morphing/operators/AddBond.hpp src/core/morphing/operators/AddBondImpl.hpp src/core/morphing/operators/AddBond.cpp include/morphing/operators/RemoveBond.hpp src/core/morphing/operators/RemoveBondImpl.hpp src/core/morphing/operators/RemoveBond.cpp include/morphing/operators/MutateAtom.hpp src/core/morphing/operators/MutateAtomImpl.hpp src/core/morphing/operators/MutateAtom.cpp include/morphing/operators/InterlayAtom.hpp src/core/morphing/operators/InterlayAtomImpl.hpp src/core/morphing/operators/InterlayAtom.cpp include/morphing/operators/ContractBond.hpp src/core/morphing/operators/ContractBondImpl.hpp src/core/morphing/operators/ContractBond.cpp include/morphing/operators/RerouteBond.hpp src/core/morphing/operators/RerouteBondImpl.hpp src/core/morphing/operators/RerouteBond.cpp)
add_library(molpher SHARED ${SOURCE_FILES})
set_target_properties(
        molpher
        PROPERTIES
            VERSION "${VERSION_NUMBER}_${BUILD_NUMBER}"
            SOVERSION "${VERSION_NUMBER_SO}"
)
message(STATUS "Libraries linked to the target molpher: ${LIBS}")
target_link_libraries(molpher ${LIBS})

## update the SWIG wrapping code if requested
set(PYTHON_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/src/python/molpher/swig_wrappers)
option(RUN_SWIG "Run SWIG" OFF)
if(RUN_SWIG)
    message(STATUS "SWIG wrapping code will be updated")
    add_custom_command(
            TARGET molpher
            PRE_BUILD
            COMMAND swig3.0 -python -I${CMAKE_SOURCE_DIR}/include/ -Wall -c++ -outdir ${PYTHON_PACKAGE_DIR} -o ${CMAKE_SOURCE_DIR}/src/swig/molpher_wrap.cpp ${CMAKE_SOURCE_DIR}/include/molpher.i
    )
else()
    message(STATUS "SWIG wrapping code will NOT be updated. Use the -DRUN_SWIG=ON to update.")
endif()

# Build the tests binary
if (CPPUNIT_LIBRARIES)
    set(TEST_SRC tests/minimal_test)
    set(TEST_HELPERS tests/helpers)
    add_executable(
            minimal_test
            EXCLUDE_FROM_ALL
            ${TEST_SRC}/MinimalTest.cpp
            ${TEST_SRC}/MinimalTest.hpp
            ${TEST_SRC}/MinimalTestRunner.cpp
            ${TEST_HELPERS}/io/stdout.hpp
            ${TEST_HELPERS}/io/stdout.cpp
            tests/helpers/mol_helpers.hpp tests/helpers/mol_helpers.cpp)
    target_include_directories(minimal_test PUBLIC ${CPPUNIT_INCLUDE_DIR})
    target_include_directories(minimal_test PUBLIC tests/helpers/)
    target_link_libraries(
            minimal_test
            molpher
            ${TEST_LIBS}
    )
endif()

# Setup installation directory (${CMAKE_SOURCE_DIR}/dist/ if not requested otherwise)
set(INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/dist/" CACHE PATH "The root of the library installation")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "CMAKE_INSTALL_PREFIX not set. Setting install path to: ${INSTALL_PREFIX}")
    SET(CMAKE_INSTALL_PREFIX
            "${INSTALL_PREFIX}" CACHE PATH "install prefix" FORCE
            )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_INSTALL_PREFIX}")

# Install Boost

## create the Boost installation target
foreach(lib ${Boost_LIBRARIES})
    list(APPEND Boost_LIBRARIES_PATTERNS "${lib}.*")
endforeach()
add_custom_command(
        OUTPUT
        ${CMAKE_INSTALL_PREFIX}/include/boost/*
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${Boost_INCLUDE_DIRS}/boost/ ${CMAKE_INSTALL_PREFIX}/include/boost
)
add_custom_command(
        OUTPUT
        ${CMAKE_INSTALL_PREFIX}/lib/libboost*.so
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib/
        COMMAND ${CMAKE_COMMAND} -E copy ${Boost_LIBRARIES} ${Boost_LIBRARIES_PATTERNS} ${CMAKE_INSTALL_PREFIX}/lib/
)

set(Boost_TARGET "boost_install")
add_custom_target(
        ${Boost_TARGET}
        DEPENDS ${CMAKE_INSTALL_PREFIX}/lib/libboost*.so ${CMAKE_INSTALL_PREFIX}/include/boost/*
        COMMENT "Installed Boost."
)

# Install RDKit

## create the RDKit installation target
foreach(lib ${RDKIT_LIBRARIES})
    list(APPEND RDKIT_LIBRARIES_PATTERNS "${lib}.*")
endforeach()
add_custom_command(
        OUTPUT
        ${CMAKE_INSTALL_PREFIX}/include/rdkit/*
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${RDKIT_INCLUDE_DIR} ${CMAKE_INSTALL_PREFIX}/include/rdkit
)
add_custom_command(
        OUTPUT
        ${CMAKE_INSTALL_PREFIX}/lib/libRDKit*.so
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib/
        COMMAND ${CMAKE_COMMAND} -E copy ${RDKIT_LIBRARIES} ${RDKIT_LIBRARIES_PATTERNS} ${CMAKE_INSTALL_PREFIX}/lib/
)

set(RDKIT_TARGET "rdkit_install")
add_custom_target(
        ${RDKIT_TARGET}
        DEPENDS ${CMAKE_INSTALL_PREFIX}/lib/libRDKit*.so ${CMAKE_INSTALL_PREFIX}/include/rdkit/*
        COMMENT "Installed RDKit."
)

# Install TBB

## create the TBB installation target
foreach(lib ${TBB_LIBRARIES})
    list(APPEND TBB_LIBRARIES_PATTERNS "${lib}.*")
endforeach()
add_custom_command(
        OUTPUT
        ${CMAKE_INSTALL_PREFIX}/include/tbb/*.h
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${TBB_INCLUDE_DIRS}/tbb/ ${CMAKE_INSTALL_PREFIX}/include/tbb
)
add_custom_command(
        OUTPUT
        ${CMAKE_INSTALL_PREFIX}/lib/libtbb*.so
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib/
        COMMAND ${CMAKE_COMMAND} -E copy ${TBB_LIBRARIES} ${TBB_LIBRARIES_PATTERNS} ${CMAKE_INSTALL_PREFIX}/lib/
)
set(TBB_TARGET "tbb_install")
add_custom_target(
        ${TBB_TARGET}
        DEPENDS ${CMAKE_INSTALL_PREFIX}/lib/libtbb*.so ${CMAKE_INSTALL_PREFIX}/include/tbb/*.h
        COMMENT "Installed TBB."
)

# Install unit tests if Cppunit available

if (CPPUNIT_LIBRARIES)
    ## create the unit test installation target
    get_filename_component(TEST_LIBS_DIR ${TEST_LIBS} DIRECTORY)
    FILE(GLOB TEST_LIBS_NAMES "${TEST_LIBS_DIR}/libcppunit*")
    add_custom_command(
            OUTPUT
            ${CMAKE_INSTALL_PREFIX}/include/cppunit/*
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${CPPUNIT_INCLUDE_DIR}/cppunit/ ${CMAKE_INSTALL_PREFIX}/include/cppunit
    )
    add_custom_command(
            OUTPUT
            ${CMAKE_INSTALL_PREFIX}/lib/libcppunit*
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib/
            COMMAND ${CMAKE_COMMAND} -E copy ${TEST_LIBS_NAMES} ${CMAKE_INSTALL_PREFIX}/lib/
    )
    add_custom_command(
            OUTPUT
            ${CMAKE_INSTALL_PREFIX}/bin/minimal_test
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/bin/
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:minimal_test> ${CMAKE_INSTALL_PREFIX}/bin/
    )
    add_custom_target(
            minimal_test_install
            DEPENDS minimal_test molpher_install ${CMAKE_INSTALL_PREFIX}/lib/libcppunit* ${CMAKE_INSTALL_PREFIX}/include/cppunit/* ${CMAKE_INSTALL_PREFIX}/bin/minimal_test
            COMMENT "Installed minimal_test."
    )
endif()

# Install Molpher

## set molpher library components to install
install(
        DIRECTORY ${CMAKE_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/molpher-lib/
        FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)
install(
        TARGETS molpher
        LIBRARY
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(
        DIRECTORY ${CMAKE_SOURCE_DIR}/res/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/molpher
)

## check for selected dependencies
set(MOLPHER_INSTALL_TARGET "molpher_install")
set(MOLPHER_INSTALL_TARGET_DEPENDS "molpher")
# TODO: make the definition of the options below more automatic
option(INSTALL_TBB "Install TBB with the library" ON)
if (INSTALL_TBB)
    message(STATUS "TBB will also be installed to ${CMAKE_INSTALL_PREFIX}. Use -DINSTALL_TBB=OFF to disable.")
    set(MOLPHER_INSTALL_TARGET_DEPENDS ${MOLPHER_INSTALL_TARGET_DEPENDS} ${TBB_TARGET})
else()
    message(STATUS "The TBB runtime dependency will NOT be installed with the library. Use -DINSTALL_TBB=ON to enable.")
endif()
option(INSTALL_Boost "Install Boost with the library" ON)
if (INSTALL_Boost)
    message(STATUS "Boost will also be installed to ${CMAKE_INSTALL_PREFIX}. Use -DINSTALL_Boost=OFF to disable.")
    set(MOLPHER_INSTALL_TARGET_DEPENDS ${MOLPHER_INSTALL_TARGET_DEPENDS} ${Boost_TARGET})
else()
    message(STATUS "The Boost runtime dependency will NOT be installed with the library. Use -DINSTALL_Boost=ON to enable.")
endif()
option(INSTALL_RDKit "Install RDKit with the library" ON)
if (INSTALL_RDKit)
    message(STATUS "RDKit will also be installed to ${CMAKE_INSTALL_PREFIX}. Use -DINSTALL_RDKit=OFF to disable.")
    set(MOLPHER_INSTALL_TARGET_DEPENDS ${MOLPHER_INSTALL_TARGET_DEPENDS} ${RDKIT_TARGET})
else()
    message(STATUS "The RDKit runtime dependency will NOT be installed with the library. Use -DINSTALL_RDKit=ON to enable.")
endif()

## add installation target
add_custom_target(${MOLPHER_INSTALL_TARGET}
        make install
        DEPENDS ${MOLPHER_INSTALL_TARGET_DEPENDS}
        COMMENT "Installing ${PROJECT_NAME} and the selected libraries...")

# Add the default installation directory to clean
set(DIST_PREFIX
        "${CMAKE_SOURCE_DIR}/dist" CACHE PATH "dist prefix" FORCE
        )
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${DIST_PREFIX}")
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_SOURCE_DIR}/build")
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_SOURCE_DIR}/src/python/molpher.egg-info")

# Build the SWIG wrappers for Python

## add target for the binaries and resource files
if(NOT PYTHON_EXECUTABLE)
    message(STATUS "PYTHON_EXECUTABLE not set. Using system defaults.")
    set(PYTHON_EXECUTABLE "python")
endif()
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;info = sys.version_info;sys.stdout.write('{0}.{1}'.format(info.major, info.minor))" OUTPUT_VARIABLE PYTHON_VERSION)
set(PYTHON_INSTALL_PREFIX
        "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages" CACHE PATH "python package install directory" FORCE
        )
set(ENV{PYTHONPATH} ${PYTHON_INSTALL_PREFIX} ENV{PYTHONPATH})
add_custom_command(
        OUTPUT
        ${PYTHON_PACKAGE_DIR}/SAScore.dat
        ${PYTHON_PACKAGE_DIR}/*.so
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/setup.py build_ext --library-dirs=${CMAKE_INSTALL_PREFIX}/lib/ --rpath=${CMAKE_INSTALL_PREFIX}/lib/
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/res/SAScore.dat ${PYTHON_PACKAGE_DIR}/SAScore.dat
        COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_INSTALL_PREFIX}
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=$ENV{PYTHONPATH} ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/setup.py install --prefix ${CMAKE_INSTALL_PREFIX}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_custom_target(molpher_install_python
        DEPENDS ${MOLPHER_INSTALL_TARGET} ${PYTHON_PACKAGE_DIR}/*.so # this works, but the library files are not cleaned for some reason...
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Built and installed SWIG wrapper module for Python...")