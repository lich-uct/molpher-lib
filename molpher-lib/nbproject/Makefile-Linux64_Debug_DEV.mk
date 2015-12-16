#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Linux64_Debug_DEV
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/molpher_API/ExplorationParameters.o \
	${OBJECTDIR}/src/molpher_API/ExplorationTree.o \
	${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot.o \
	${OBJECTDIR}/src/molpher_API/MolpherMol.o \
	${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback.o \
	${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback.o \
	${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper.o \
	${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper.o \
	${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper.o \
	${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper.o \
	${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper.o \
	${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper.o \
	${OBJECTDIR}/src/molpher_API/operations/TraverseOper.o \
	${OBJECTDIR}/src/molpher_API/operations/TreeOperation.o \
	${OBJECTDIR}/src/morphing_functions.o \
	${OBJECTDIR}/src/path_finders/BasicPathFinder.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -DNOT_NETBEANS -Wno-deprecated -Wno-write-strings -Wno-attributes -Wno-strict-aliasing -fpermissive -fPIC
CXXFLAGS=-m64 -DNOT_NETBEANS -Wno-deprecated -Wno-write-strings -Wno-attributes -Wno-strict-aliasing -fpermissive -fPIC

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L../backend/dist/Linux64_Release_library/GNU_4_6-Linux-x86 ../backend/dist/Linux64_Debug_library/GNU_4_6-Linux-x86/libbackend.a ../dependencies/rdkit/lib/libDistGeomHelpers_static.a ../dependencies/rdkit/lib/libMolAlign_static.a ../dependencies/rdkit/lib/libAlignment_static.a ../dependencies/rdkit/lib/libFragCatalog_static.a ../dependencies/rdkit/lib/libMolCatalog_static.a ../dependencies/rdkit/lib/libCatalogs_static.a ../dependencies/rdkit/lib/libChemReactions_static.a ../dependencies/rdkit/lib/libDepictor_static.a ../dependencies/rdkit/lib/libDescriptors_static.a ../dependencies/rdkit/lib/libDistGeometry_static.a ../dependencies/rdkit/lib/libFileParsers_static.a ../dependencies/rdkit/lib/libFingerprints_static.a ../dependencies/rdkit/lib/libSubgraphs_static.a ../dependencies/rdkit/lib/libForceFieldHelpers_static.a ../dependencies/rdkit/lib/libChemTransforms_static.a ../dependencies/rdkit/lib/libMolChemicalFeatures_static.a ../dependencies/rdkit/lib/libSubstructMatch_static.a ../dependencies/rdkit/lib/libSmilesParse_static.a ../dependencies/rdkit/lib/libPartialCharges_static.a ../dependencies/rdkit/lib/libShapeHelpers_static.a ../dependencies/rdkit/lib/libMolTransforms_static.a ../dependencies/rdkit/lib/libSLNParse_static.a ../dependencies/rdkit/lib/libGraphMol_static.a ../dependencies/rdkit/lib/libForceField_static.a ../dependencies/rdkit/lib/libRDGeometryLib_static.a ../dependencies/rdkit/lib/libDataStructs_static.a ../dependencies/rdkit/lib/libOptimizer_static.a ../dependencies/rdkit/lib/libEigenSolvers_static.a ../dependencies/rdkit/lib/libChemicalFeatures_static.a ../dependencies/rdkit/lib/libSimDivPickers_static.a ../dependencies/rdkit/lib/libRDGeneral_static.a ../dependencies/rdkit/lib/libhc_static.a ../dependencies/rcf/bin/libRcfLib.a ../dependencies/boost/stage/lib/libboost_date_time.a ../dependencies/boost/stage/lib/libboost_filesystem.a ../dependencies/boost/stage/lib/libboost_program_options.a ../dependencies/boost/stage/lib/libboost_regex.a ../dependencies/boost/stage/lib/libboost_wserialization.a ../dependencies/boost/stage/lib/libboost_serialization.a ../dependencies/boost/stage/lib/libboost_signals.a ../dependencies/boost/stage/lib/libboost_thread.a ../dependencies/boost/stage/lib/libboost_system.a ../dependencies/zlib/libz.a ../dependencies/tbb/lib/intel64/gcc4.4/libtbb.so.2 ../dependencies/tbb/lib/intel64/gcc4.4/libtbbmalloc_proxy.so.2 ../dependencies/tbb/lib/intel64/gcc4.4/libtbbmalloc.so.2

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../backend/dist/Linux64_Debug_library/GNU_4_6-Linux-x86/libbackend.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libDistGeomHelpers_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libMolAlign_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libAlignment_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libFragCatalog_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libMolCatalog_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libCatalogs_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libChemReactions_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libDepictor_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libDescriptors_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libDistGeometry_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libFileParsers_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libFingerprints_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libSubgraphs_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libForceFieldHelpers_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libChemTransforms_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libMolChemicalFeatures_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libSubstructMatch_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libSmilesParse_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libPartialCharges_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libShapeHelpers_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libMolTransforms_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libSLNParse_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libGraphMol_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libForceField_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libRDGeometryLib_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libDataStructs_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libOptimizer_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libEigenSolvers_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libChemicalFeatures_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libSimDivPickers_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libRDGeneral_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rdkit/lib/libhc_static.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/rcf/bin/libRcfLib.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_date_time.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_filesystem.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_program_options.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_regex.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_wserialization.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_serialization.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_signals.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_thread.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/boost/stage/lib/libboost_system.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/zlib/libz.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/tbb/lib/intel64/gcc4.4/libtbb.so.2

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/tbb/lib/intel64/gcc4.4/libtbbmalloc_proxy.so.2

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ../dependencies/tbb/lib/intel64/gcc4.4/libtbbmalloc.so.2

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -lpthread -Wl,-rpath,'$$ORIGIN/' -shared -fPIC

${OBJECTDIR}/src/molpher_API/ExplorationParameters.o: src/molpher_API/ExplorationParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/ExplorationParameters.o src/molpher_API/ExplorationParameters.cpp

${OBJECTDIR}/src/molpher_API/ExplorationTree.o: src/molpher_API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/ExplorationTree.o src/molpher_API/ExplorationTree.cpp

${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot.o: src/molpher_API/ExplorationTreeSnapshot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot.o src/molpher_API/ExplorationTreeSnapshot.cpp

${OBJECTDIR}/src/molpher_API/MolpherMol.o: src/molpher_API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/MolpherMol.o src/molpher_API/MolpherMol.cpp

${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback.o: src/molpher_API/callbacks/EraseSubtreeCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback.o src/molpher_API/callbacks/EraseSubtreeCallback.cpp

${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback.o: src/molpher_API/callbacks/TraverseCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback.o src/molpher_API/callbacks/TraverseCallback.cpp

${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper.o: src/molpher_API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper.o src/molpher_API/operations/ExtendTreeOper.cpp

${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper.o: src/molpher_API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper.o src/molpher_API/operations/FilterMorphsOper.cpp

${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper.o: src/molpher_API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper.o src/molpher_API/operations/FindLeavesOper.cpp

${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper.o: src/molpher_API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper.o src/molpher_API/operations/GenerateMorphsOper.cpp

${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper.o: src/molpher_API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper.o src/molpher_API/operations/PruneTreeOper.cpp

${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper.o: src/molpher_API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper.o src/molpher_API/operations/SortMorphsOper.cpp

${OBJECTDIR}/src/molpher_API/operations/TraverseOper.o: src/molpher_API/operations/TraverseOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/TraverseOper.o src/molpher_API/operations/TraverseOper.cpp

${OBJECTDIR}/src/molpher_API/operations/TreeOperation.o: src/molpher_API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/TreeOperation.o src/molpher_API/operations/TreeOperation.cpp

${OBJECTDIR}/src/morphing_functions.o: src/morphing_functions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/morphing_functions.o src/morphing_functions.cpp

${OBJECTDIR}/src/path_finders/BasicPathFinder.o: src/path_finders/BasicPathFinder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/path_finders
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/path_finders/BasicPathFinder.o src/path_finders/BasicPathFinder.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/APITests.o ${TESTDIR}/tests/testrunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}  -pthread -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} -L/usr/lib/python3.4/config-3.4m-x86_64-linux-gnu -Wl,-rpath,../dependencies/tbb/lib/intel64/gcc4.4 `cppunit-config --libs` /usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4.so   


${TESTDIR}/tests/APITests.o: tests/APITests.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/APITests.o tests/APITests.cpp


${TESTDIR}/tests/testrunner.o: tests/testrunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/testrunner.o tests/testrunner.cpp


${OBJECTDIR}/src/molpher_API/ExplorationParameters_nomain.o: ${OBJECTDIR}/src/molpher_API/ExplorationParameters.o src/molpher_API/ExplorationParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/ExplorationParameters.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/ExplorationParameters_nomain.o src/molpher_API/ExplorationParameters.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/ExplorationParameters.o ${OBJECTDIR}/src/molpher_API/ExplorationParameters_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/ExplorationTree_nomain.o: ${OBJECTDIR}/src/molpher_API/ExplorationTree.o src/molpher_API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/ExplorationTree.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/ExplorationTree_nomain.o src/molpher_API/ExplorationTree.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/ExplorationTree.o ${OBJECTDIR}/src/molpher_API/ExplorationTree_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot_nomain.o: ${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot.o src/molpher_API/ExplorationTreeSnapshot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot_nomain.o src/molpher_API/ExplorationTreeSnapshot.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot.o ${OBJECTDIR}/src/molpher_API/ExplorationTreeSnapshot_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/MolpherMol_nomain.o: ${OBJECTDIR}/src/molpher_API/MolpherMol.o src/molpher_API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/MolpherMol.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/MolpherMol_nomain.o src/molpher_API/MolpherMol.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/MolpherMol.o ${OBJECTDIR}/src/molpher_API/MolpherMol_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback_nomain.o: ${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback.o src/molpher_API/callbacks/EraseSubtreeCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback_nomain.o src/molpher_API/callbacks/EraseSubtreeCallback.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback.o ${OBJECTDIR}/src/molpher_API/callbacks/EraseSubtreeCallback_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback_nomain.o: ${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback.o src/molpher_API/callbacks/TraverseCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback_nomain.o src/molpher_API/callbacks/TraverseCallback.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback.o ${OBJECTDIR}/src/molpher_API/callbacks/TraverseCallback_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper.o src/molpher_API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper_nomain.o src/molpher_API/operations/ExtendTreeOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper.o ${OBJECTDIR}/src/molpher_API/operations/ExtendTreeOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper.o src/molpher_API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper_nomain.o src/molpher_API/operations/FilterMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper.o ${OBJECTDIR}/src/molpher_API/operations/FilterMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper.o src/molpher_API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper_nomain.o src/molpher_API/operations/FindLeavesOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper.o ${OBJECTDIR}/src/molpher_API/operations/FindLeavesOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper.o src/molpher_API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper_nomain.o src/molpher_API/operations/GenerateMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper.o ${OBJECTDIR}/src/molpher_API/operations/GenerateMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper.o src/molpher_API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper_nomain.o src/molpher_API/operations/PruneTreeOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper.o ${OBJECTDIR}/src/molpher_API/operations/PruneTreeOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper.o src/molpher_API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper_nomain.o src/molpher_API/operations/SortMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper.o ${OBJECTDIR}/src/molpher_API/operations/SortMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/TraverseOper_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/TraverseOper.o src/molpher_API/operations/TraverseOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/TraverseOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/TraverseOper_nomain.o src/molpher_API/operations/TraverseOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/TraverseOper.o ${OBJECTDIR}/src/molpher_API/operations/TraverseOper_nomain.o;\
	fi

${OBJECTDIR}/src/molpher_API/operations/TreeOperation_nomain.o: ${OBJECTDIR}/src/molpher_API/operations/TreeOperation.o src/molpher_API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/molpher_API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/molpher_API/operations/TreeOperation.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/molpher_API/operations/TreeOperation_nomain.o src/molpher_API/operations/TreeOperation.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/molpher_API/operations/TreeOperation.o ${OBJECTDIR}/src/molpher_API/operations/TreeOperation_nomain.o;\
	fi

${OBJECTDIR}/src/morphing_functions_nomain.o: ${OBJECTDIR}/src/morphing_functions.o src/morphing_functions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/morphing_functions.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/morphing_functions_nomain.o src/morphing_functions.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/morphing_functions.o ${OBJECTDIR}/src/morphing_functions_nomain.o;\
	fi

${OBJECTDIR}/src/path_finders/BasicPathFinder_nomain.o: ${OBJECTDIR}/src/path_finders/BasicPathFinder.o src/path_finders/BasicPathFinder.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/path_finders
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/path_finders/BasicPathFinder.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./include -I../common -I../backend/ -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/path_finders/BasicPathFinder_nomain.o src/path_finders/BasicPathFinder.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/path_finders/BasicPathFinder.o ${OBJECTDIR}/src/path_finders/BasicPathFinder_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
