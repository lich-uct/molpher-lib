cmake_minimum_required(VERSION 3.5)
project(molpher_lib)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set some general tool flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath,$ORIGIN/")

# search paths for the dependencies
set(DEPS_DIR "deps")
set(
        CMAKE_INCLUDE_PATH

        ${CMAKE_INCLUDE_PATH}
        ${DEPS_DIR}/boost/
        ${DEPS_DIR}/rdkit/Code/
        ${DEPS_DIR}/tbb/include/
        ENV{CPATH}
)
set(
        CMAKE_LIBRARY_PATH

        ${CMAKE_LIBRARY_PATH}
        ${DEPS_DIR}/boost/stage/lib/
        ${DEPS_DIR}/rdkit/lib/
        ${DEPS_DIR}/tbb/lib/intel64/gcc4.4/
        ENV{LD_LIBRARY_PATH}
)

# include project directories
include_directories(
        include/
        src/
)

# Find libraries

## Threads
find_package(Threads REQUIRED)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "Threads library not found")
endif()
include_directories(${Threads_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Threads_LIBRARIES})
add_definitions(${Threads_DEFINITIONS})

## Boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost library not found")
endif()
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})
add_definitions(${Boost_DEFINITIONS})
add_definitions("-DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB")

## RDKit
set(RDKIT_LINK_STATIC ON)
find_package(RDKit REQUIRED)
if(NOT RDKIT_FOUND)
    message(FATAL_ERROR "RDKit library not found")
endif()
include_directories(${RDKIT_INCLUDE_DIR})
set(LIBS ${LIBS} ${RDKIT_LIBRARIES})
add_definitions(${RDKIT_DEFINITIONS})

## TBB
find_package(TBB REQUIRED)
if(NOT TBB_FOUND)
    message(FATAL_ERROR "TBB library not found")
endif()
include_directories(${TBB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${TBB_LIBRARIES})
add_definitions(${TBB_DEFINITIONS})

# print the variables set so far
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(WARNING "${_variableName}=${${_variableName}}")
#endforeach()

# Target definition
# Molpher-lib source files
set(SOURCE_FILES
        include/data_structs/ExplorationData.hpp
        include/data_structs/ExplorationTree.hpp
        include/data_structs/MolpherMol.hpp
        include/operations/callbacks/TraverseCallback.hpp
        include/operations/ExtendTreeOper.hpp
        include/operations/FilterMorphsOper.hpp
        include/operations/FindLeavesOper.hpp
        include/operations/GenerateMorphsOper.hpp
        include/operations/PruneTreeOper.hpp
        include/operations/SortMorphsOper.hpp
        include/operations/TraverseOper.hpp
        include/operations/TreeOperation.hpp
        include/selectors/chemoper_selectors.h
        include/selectors/fingerprint_selectors.h
        include/selectors/simcoeff_selectors.h
        include/SAScore_data_loader.hpp
        src/core/API/operations/callbacks/TraverseCallback.cpp
        src/core/API/operations/callbacks/TraverseCallbackImpl.hpp
        src/core/API/operations/ExtendTreeOper.cpp
        src/core/API/operations/ExtendTreeOperImpl.hpp
        src/core/API/operations/FilterMorphsOper.cpp
        src/core/API/operations/FilterMorphsOperImpl.hpp
        src/core/API/operations/FindLeavesOper.cpp
        src/core/API/operations/FindLeavesOperImpl.hpp
        src/core/API/operations/GenerateMorphsOper.cpp
        src/core/API/operations/GenerateMorphsOperImpl.hpp
        src/core/API/operations/PruneTreeOper.cpp
        src/core/API/operations/PruneTreeOperImpl.hpp
        src/core/API/operations/SortMorphsOper.cpp
        src/core/API/operations/SortMorphsOperImpl.hpp
        src/core/API/operations/TraverseOper.cpp
        src/core/API/operations/TraverseOperImpl.hpp
        src/core/API/operations/TreeOperation.cpp
        src/core/API/operations/TreeOperationImpl.hpp
        src/core/API/ExplorationData.cpp
        src/core/API/ExplorationTree.cpp
        src/core/API/ExplorationTreeImpl.h
        src/core/API/MolpherMol.cpp
        src/core/API/MolpherMolImpl.hpp
        src/core/chem/fingerprintStrategy/AtomPairsFngpr.cpp
        src/core/chem/fingerprintStrategy/AtomPairsFngpr.hpp
        src/core/chem/fingerprintStrategy/FingerprintStrategy.cpp
        src/core/chem/fingerprintStrategy/FingerprintStrategy.h
        src/core/chem/fingerprintStrategy/MorganFngpr.cpp
        src/core/chem/fingerprintStrategy/MorganFngpr.hpp
        src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp
        src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.hpp
        src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp
        src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.hpp
        src/core/chem/fingerprintStrategy/TopolSingleFngpr.cpp
        src/core/chem/fingerprintStrategy/TopolSingleFngpr.hpp
        src/core/chem/fingerprintStrategy/TopolTorsFngpr.cpp
        src/core/chem/fingerprintStrategy/TopolTorsFngpr.hpp
        src/core/chem/fingerprintStrategy/VectorFpFngpr.cpp
        src/core/chem/fingerprintStrategy/VectorFpFngpr.hpp
        src/core/chem/morphing/CalculateDistances.hpp
        src/core/chem/morphing/CalculateMorphs.hpp
        src/core/chem/morphing/Morphing.cpp
        src/core/chem/morphing/Morphing.hpp
        src/core/chem/morphing/MorphingData.cpp
        src/core/chem/morphing/MorphingData.h
        src/core/chem/morphing/MorphingFtors.cpp
        src/core/chem/morphing/MorphingFtors.hpp
        src/core/chem/morphing/ReturnResults.hpp
        src/core/chem/morphingStrategy/MorphingStrategy.h
        src/core/chem/morphingStrategy/OpAddAtom.cpp
        src/core/chem/morphingStrategy/OpAddAtom.hpp
        src/core/chem/morphingStrategy/OpAddBond.cpp
        src/core/chem/morphingStrategy/OpAddBond.hpp
        src/core/chem/morphingStrategy/OpBondContraction.cpp
        src/core/chem/morphingStrategy/OpBondContraction.hpp
        src/core/chem/morphingStrategy/OpBondReroute.cpp
        src/core/chem/morphingStrategy/OpBondReroute.hpp
        src/core/chem/morphingStrategy/OpInterlayAtom.cpp
        src/core/chem/morphingStrategy/OpInterlayAtom.hpp
        src/core/chem/morphingStrategy/OpMutateAtom.cpp
        src/core/chem/morphingStrategy/OpMutateAtom.hpp
        src/core/chem/morphingStrategy/OpRemoveAtom.cpp
        src/core/chem/morphingStrategy/OpRemoveAtom.hpp
        src/core/chem/morphingStrategy/OpRemoveBond.cpp
        src/core/chem/morphingStrategy/OpRemoveBond.hpp
        src/core/chem/simCoefStrategy/AllBitSimCoef.cpp
        src/core/chem/simCoefStrategy/AllBitSimCoef.hpp
        src/core/chem/simCoefStrategy/AsymmetricSimCoef.cpp
        src/core/chem/simCoefStrategy/AsymmetricSimCoef.hpp
        src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp
        src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.hpp
        src/core/chem/simCoefStrategy/CosineSimCoef.cpp
        src/core/chem/simCoefStrategy/CosineSimCoef.hpp
        src/core/chem/simCoefStrategy/DiceSimCoef.cpp
        src/core/chem/simCoefStrategy/DiceSimCoef.hpp
        src/core/chem/simCoefStrategy/KulczynskiSimCoef.cpp
        src/core/chem/simCoefStrategy/KulczynskiSimCoef.hpp
        src/core/chem/simCoefStrategy/McConnaugheySimCoef.cpp
        src/core/chem/simCoefStrategy/McConnaugheySimCoef.hpp
        src/core/chem/simCoefStrategy/OnBitSimCoef.cpp
        src/core/chem/simCoefStrategy/OnBitSimCoef.hpp
        src/core/chem/simCoefStrategy/RusselSimCoef.cpp
        src/core/chem/simCoefStrategy/RusselSimCoef.hpp
        src/core/chem/simCoefStrategy/SimCoefStrategy.h
        src/core/chem/simCoefStrategy/SokalSimCoef.cpp
        src/core/chem/simCoefStrategy/SokalSimCoef.hpp
        src/core/chem/simCoefStrategy/TanimotoSimCoef.cpp
        src/core/chem/simCoefStrategy/TanimotoSimCoef.hpp
        src/core/chem/simCoefStrategy/TverskySimCoef.cpp
        src/core/chem/simCoefStrategy/TverskySimCoef.hpp
        src/core/chem/ChemicalAuxiliary.cpp
        src/core/chem/ChemicalAuxiliary.h
        src/core/chem/SimCoefCalculator.cpp
        src/core/chem/SimCoefCalculator.hpp
        src/core/data_structs/ExplorationDataImpl.hpp
        src/core/data_structs/MolpherAtom.h
        src/core/data_structs/MolpherMolData.hpp
        src/core/data_structs/MolpherParam.h
        src/core/misc/selectors/chemoper_selectors.cpp
        src/core/misc/selectors/fingerprint_selectors.cpp
        src/core/misc/selectors/simcoeff_selectors.cpp
        src/core/misc/global_types.h
        src/core/misc/inout.cpp
        src/core/misc/inout.h
        src/core/misc/iteration_serializer.cpp
        src/core/misc/iteration_serializer.hpp
        src/core/misc/SAScore.cpp
        src/core/misc/SAScore.h
        src/core/misc/SAScore_data_loader.cpp
        src/core/misc/SynchRand.cpp
        src/core/misc/SynchRand.h
        )
add_library(molpher SHARED ${SOURCE_FILES})

# Linked libraries
target_link_libraries(molpher ${LIBS})

# SWIG
#set(PYTHON_PACKAGE_DIR ${CMAKE_SOURCE_DIR}/src/python/molpher/swig_wrappers)
#option(RUN_SWIG "Run SWIG" ON)
#if(RUN_SWIG)
#    add_custom_command(
#            TARGET molpher
#            PRE_BUILD
#            COMMAND swig3.0 -python -Iinclude/ -Wall -c++ -outdir ${CMAKE_SOURCE_DIR}/src/python/molpher/swig_wrappers -o ${CMAKE_SOURCE_DIR}/src/swig/molpher_wrap.cpp ${CMAKE_SOURCE_DIR}/include/molpher.i
#    )
#endif()
#
## copy the built library to the dist folder
#add_custom_command(
#        TARGET molpher
#        POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E make_directory dist/lib/
#        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/* dist/lib/
#)
#
## copy resource files
#configure_file(${CMAKE_SOURCE_DIR}/res/SAScore.dat ${PYTHON_PACKAGE_DIR}/SAScore.dat COPYONLY)