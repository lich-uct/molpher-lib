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
CND_CONF=Debug_DEV
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/core/API/ExplorationData.o \
	${OBJECTDIR}/src/core/API/ExplorationTree.o \
	${OBJECTDIR}/src/core/API/MolpherMol.o \
	${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o \
	${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o \
	${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o \
	${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o \
	${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o \
	${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o \
	${OBJECTDIR}/src/core/API/operations/TreeOperation.o \
	${OBJECTDIR}/src/core/chem/ChemicalAuxiliary.o \
	${OBJECTDIR}/src/core/chem/SimCoefCalculator.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr.o \
	${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr.o \
	${OBJECTDIR}/src/core/chem/morphing/Morphing.o \
	${OBJECTDIR}/src/core/chem/morphing/MorphingData.o \
	${OBJECTDIR}/src/core/chem/morphing/MorphingFtors.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom.o \
	${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef.o \
	${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef.o \
	${OBJECTDIR}/src/core/misc/SAScore.o \
	${OBJECTDIR}/src/core/misc/SAScore_data_loader.o \
	${OBJECTDIR}/src/core/misc/SynchRand.o \
	${OBJECTDIR}/src/core/misc/inout.o \
	${OBJECTDIR}/src/core/misc/iteration_serializer.o \
	${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors.o \
	${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o \
	${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f2

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=deps/rdkit/lib/libDistGeomHelpers_static.a deps/rdkit/lib/libMolAlign_static.a deps/rdkit/lib/libAlignment_static.a deps/rdkit/lib/libFragCatalog_static.a deps/rdkit/lib/libMolCatalog_static.a deps/rdkit/lib/libCatalogs_static.a deps/rdkit/lib/libChemReactions_static.a deps/rdkit/lib/libDepictor_static.a deps/rdkit/lib/libDescriptors_static.a deps/rdkit/lib/libDistGeometry_static.a deps/rdkit/lib/libFileParsers_static.a deps/rdkit/lib/libFingerprints_static.a deps/rdkit/lib/libSubgraphs_static.a deps/rdkit/lib/libForceFieldHelpers_static.a deps/rdkit/lib/libChemTransforms_static.a deps/rdkit/lib/libMolChemicalFeatures_static.a deps/rdkit/lib/libSubstructMatch_static.a deps/rdkit/lib/libSmilesParse_static.a deps/rdkit/lib/libPartialCharges_static.a deps/rdkit/lib/libShapeHelpers_static.a deps/rdkit/lib/libMolTransforms_static.a deps/rdkit/lib/libSLNParse_static.a deps/rdkit/lib/libGraphMol_static.a deps/rdkit/lib/libForceField_static.a deps/rdkit/lib/libRDGeometryLib_static.a deps/rdkit/lib/libDataStructs_static.a deps/rdkit/lib/libOptimizer_static.a deps/rdkit/lib/libEigenSolvers_static.a deps/rdkit/lib/libChemicalFeatures_static.a deps/rdkit/lib/libSimDivPickers_static.a deps/rdkit/lib/libRDGeneral_static.a deps/rdkit/lib/libhc_static.a deps/boost/stage/lib/libboost_date_time.a deps/boost/stage/lib/libboost_filesystem.a deps/boost/stage/lib/libboost_program_options.a deps/boost/stage/lib/libboost_regex.a deps/boost/stage/lib/libboost_wserialization.a deps/boost/stage/lib/libboost_serialization.a deps/boost/stage/lib/libboost_signals.a deps/boost/stage/lib/libboost_thread.a deps/boost/stage/lib/libboost_system.a deps/tbb/lib/intel64/gcc4.4/libtbbmalloc_proxy_debug.so.2 deps/tbb/lib/intel64/gcc4.4/libtbb_debug.so.2 deps/tbb/lib/intel64/gcc4.4/libtbbmalloc_debug.so.2

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk dist/lib/libmolpher.so

dist/lib/libmolpher.so: deps/rdkit/lib/libDistGeomHelpers_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libMolAlign_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libAlignment_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libFragCatalog_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libMolCatalog_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libCatalogs_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libChemReactions_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libDepictor_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libDescriptors_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libDistGeometry_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libFileParsers_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libFingerprints_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libSubgraphs_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libForceFieldHelpers_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libChemTransforms_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libMolChemicalFeatures_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libSubstructMatch_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libSmilesParse_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libPartialCharges_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libShapeHelpers_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libMolTransforms_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libSLNParse_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libGraphMol_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libForceField_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libRDGeometryLib_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libDataStructs_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libOptimizer_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libEigenSolvers_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libChemicalFeatures_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libSimDivPickers_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libRDGeneral_static.a

dist/lib/libmolpher.so: deps/rdkit/lib/libhc_static.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_date_time.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_filesystem.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_program_options.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_regex.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_wserialization.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_serialization.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_signals.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_thread.a

dist/lib/libmolpher.so: deps/boost/stage/lib/libboost_system.a

dist/lib/libmolpher.so: deps/tbb/lib/intel64/gcc4.4/libtbbmalloc_proxy_debug.so.2

dist/lib/libmolpher.so: deps/tbb/lib/intel64/gcc4.4/libtbb_debug.so.2

dist/lib/libmolpher.so: deps/tbb/lib/intel64/gcc4.4/libtbbmalloc_debug.so.2

dist/lib/libmolpher.so: ${OBJECTFILES}
	${MKDIR} -p dist/lib
	${LINK.cc} -o dist/lib/libmolpher.so ${OBJECTFILES} ${LDLIBSOPTIONS} -lpthread -Wl,-rpath,'$$ORIGIN/' -shared -fPIC

${OBJECTDIR}/src/core/API/ExplorationData.o: src/core/API/ExplorationData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationData.o src/core/API/ExplorationData.cpp

${OBJECTDIR}/src/core/API/ExplorationTree.o: src/core/API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationTree.o src/core/API/ExplorationTree.cpp

${OBJECTDIR}/src/core/API/MolpherMol.o: src/core/API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/MolpherMol.o src/core/API/MolpherMol.cpp

${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o: src/core/API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o src/core/API/operations/ExtendTreeOper.cpp

${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o: src/core/API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o src/core/API/operations/FilterMorphsOper.cpp

${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o: src/core/API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o src/core/API/operations/FindLeavesOper.cpp

${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o: src/core/API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o src/core/API/operations/GenerateMorphsOper.cpp

${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o: src/core/API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o src/core/API/operations/PruneTreeOper.cpp

${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o: src/core/API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o src/core/API/operations/SortMorphsOper.cpp

${OBJECTDIR}/src/core/API/operations/TreeOperation.o: src/core/API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/TreeOperation.o src/core/API/operations/TreeOperation.cpp

${OBJECTDIR}/src/core/chem/ChemicalAuxiliary.o: src/core/chem/ChemicalAuxiliary.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/ChemicalAuxiliary.o src/core/chem/ChemicalAuxiliary.cpp

${OBJECTDIR}/src/core/chem/SimCoefCalculator.o: src/core/chem/SimCoefCalculator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/SimCoefCalculator.o src/core/chem/SimCoefCalculator.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr.o: src/core/chem/fingerprintStrategy/AtomPairsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr.o src/core/chem/fingerprintStrategy/AtomPairsFngpr.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy.o: src/core/chem/fingerprintStrategy/FingerprintStrategy.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy.o src/core/chem/fingerprintStrategy/FingerprintStrategy.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr.o: src/core/chem/fingerprintStrategy/MorganFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr.o src/core/chem/fingerprintStrategy/MorganFngpr.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.o: src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.o src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.o: src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.o src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr.o: src/core/chem/fingerprintStrategy/TopolSingleFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr.o src/core/chem/fingerprintStrategy/TopolSingleFngpr.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr.o: src/core/chem/fingerprintStrategy/TopolTorsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr.o src/core/chem/fingerprintStrategy/TopolTorsFngpr.cpp

${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr.o: src/core/chem/fingerprintStrategy/VectorFpFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr.o src/core/chem/fingerprintStrategy/VectorFpFngpr.cpp

${OBJECTDIR}/src/core/chem/morphing/Morphing.o: src/core/chem/morphing/Morphing.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphing/Morphing.o src/core/chem/morphing/Morphing.cpp

${OBJECTDIR}/src/core/chem/morphing/MorphingData.o: src/core/chem/morphing/MorphingData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphing/MorphingData.o src/core/chem/morphing/MorphingData.cpp

${OBJECTDIR}/src/core/chem/morphing/MorphingFtors.o: src/core/chem/morphing/MorphingFtors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphing/MorphingFtors.o src/core/chem/morphing/MorphingFtors.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom.o: src/core/chem/morphingStrategy/OpAddAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom.o src/core/chem/morphingStrategy/OpAddAtom.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond.o: src/core/chem/morphingStrategy/OpAddBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond.o src/core/chem/morphingStrategy/OpAddBond.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction.o: src/core/chem/morphingStrategy/OpBondContraction.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction.o src/core/chem/morphingStrategy/OpBondContraction.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute.o: src/core/chem/morphingStrategy/OpBondReroute.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute.o src/core/chem/morphingStrategy/OpBondReroute.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom.o: src/core/chem/morphingStrategy/OpInterlayAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom.o src/core/chem/morphingStrategy/OpInterlayAtom.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom.o: src/core/chem/morphingStrategy/OpMutateAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom.o src/core/chem/morphingStrategy/OpMutateAtom.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom.o: src/core/chem/morphingStrategy/OpRemoveAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom.o src/core/chem/morphingStrategy/OpRemoveAtom.cpp

${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond.o: src/core/chem/morphingStrategy/OpRemoveBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond.o src/core/chem/morphingStrategy/OpRemoveBond.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef.o: src/core/chem/simCoefStrategy/AllBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef.o src/core/chem/simCoefStrategy/AllBitSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef.o: src/core/chem/simCoefStrategy/AsymmetricSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef.o src/core/chem/simCoefStrategy/AsymmetricSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.o: src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.o src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef.o: src/core/chem/simCoefStrategy/CosineSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef.o src/core/chem/simCoefStrategy/CosineSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef.o: src/core/chem/simCoefStrategy/DiceSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef.o src/core/chem/simCoefStrategy/DiceSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef.o: src/core/chem/simCoefStrategy/KulczynskiSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef.o src/core/chem/simCoefStrategy/KulczynskiSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef.o: src/core/chem/simCoefStrategy/McConnaugheySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef.o src/core/chem/simCoefStrategy/McConnaugheySimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef.o: src/core/chem/simCoefStrategy/OnBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef.o src/core/chem/simCoefStrategy/OnBitSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef.o: src/core/chem/simCoefStrategy/RusselSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef.o src/core/chem/simCoefStrategy/RusselSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef.o: src/core/chem/simCoefStrategy/SokalSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef.o src/core/chem/simCoefStrategy/SokalSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef.o: src/core/chem/simCoefStrategy/TanimotoSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef.o src/core/chem/simCoefStrategy/TanimotoSimCoef.cpp

${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef.o: src/core/chem/simCoefStrategy/TverskySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef.o src/core/chem/simCoefStrategy/TverskySimCoef.cpp

${OBJECTDIR}/src/core/misc/SAScore.o: src/core/misc/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SAScore.o src/core/misc/SAScore.cpp

${OBJECTDIR}/src/core/misc/SAScore_data_loader.o: src/core/misc/SAScore_data_loader.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SAScore_data_loader.o src/core/misc/SAScore_data_loader.cpp

${OBJECTDIR}/src/core/misc/SynchRand.o: src/core/misc/SynchRand.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SynchRand.o src/core/misc/SynchRand.cpp

${OBJECTDIR}/src/core/misc/inout.o: src/core/misc/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/inout.o src/core/misc/inout.cpp

${OBJECTDIR}/src/core/misc/iteration_serializer.o: src/core/misc/iteration_serializer.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/iteration_serializer.o src/core/misc/iteration_serializer.cpp

${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors.o: src/core/misc/selectors/chemoper_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors.o src/core/misc/selectors/chemoper_selectors.cpp

${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o: src/core/misc/selectors/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o src/core/misc/selectors/fingerprint_selectors.cpp

${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o: src/core/misc/selectors/simcoeff_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o src/core/misc/selectors/simcoeff_selectors.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/minimal_test/MinimalTest.o ${TESTDIR}/tests/minimal_test/MinimalTestRunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}  -pthread -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} -Wl,-rpath,dist/lib `cppunit-config --libs`   


${TESTDIR}/tests/minimal_test/MinimalTest.o: tests/minimal_test/MinimalTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests/minimal_test
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Iinclude/ -std=c++11 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/minimal_test/MinimalTest.o tests/minimal_test/MinimalTest.cpp


${TESTDIR}/tests/minimal_test/MinimalTestRunner.o: tests/minimal_test/MinimalTestRunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests/minimal_test
	${RM} "$@.d"
	$(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Iinclude/ -std=c++11 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/minimal_test/MinimalTestRunner.o tests/minimal_test/MinimalTestRunner.cpp


${OBJECTDIR}/src/core/API/ExplorationData_nomain.o: ${OBJECTDIR}/src/core/API/ExplorationData.o src/core/API/ExplorationData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/ExplorationData.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationData_nomain.o src/core/API/ExplorationData.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/ExplorationData.o ${OBJECTDIR}/src/core/API/ExplorationData_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/ExplorationTree_nomain.o: ${OBJECTDIR}/src/core/API/ExplorationTree.o src/core/API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/ExplorationTree.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationTree_nomain.o src/core/API/ExplorationTree.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/ExplorationTree.o ${OBJECTDIR}/src/core/API/ExplorationTree_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/MolpherMol_nomain.o: ${OBJECTDIR}/src/core/API/MolpherMol.o src/core/API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/MolpherMol.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/MolpherMol_nomain.o src/core/API/MolpherMol.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/MolpherMol.o ${OBJECTDIR}/src/core/API/MolpherMol_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/ExtendTreeOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o src/core/API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper_nomain.o src/core/API/operations/ExtendTreeOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/FilterMorphsOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o src/core/API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper_nomain.o src/core/API/operations/FilterMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/FindLeavesOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o src/core/API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FindLeavesOper_nomain.o src/core/API/operations/FindLeavesOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o ${OBJECTDIR}/src/core/API/operations/FindLeavesOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o src/core/API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper_nomain.o src/core/API/operations/GenerateMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/PruneTreeOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o src/core/API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/PruneTreeOper_nomain.o src/core/API/operations/PruneTreeOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o ${OBJECTDIR}/src/core/API/operations/PruneTreeOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/SortMorphsOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o src/core/API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/SortMorphsOper_nomain.o src/core/API/operations/SortMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o ${OBJECTDIR}/src/core/API/operations/SortMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/TreeOperation_nomain.o: ${OBJECTDIR}/src/core/API/operations/TreeOperation.o src/core/API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/TreeOperation.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/TreeOperation_nomain.o src/core/API/operations/TreeOperation.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/TreeOperation.o ${OBJECTDIR}/src/core/API/operations/TreeOperation_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/ChemicalAuxiliary_nomain.o: ${OBJECTDIR}/src/core/chem/ChemicalAuxiliary.o src/core/chem/ChemicalAuxiliary.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/ChemicalAuxiliary.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/ChemicalAuxiliary_nomain.o src/core/chem/ChemicalAuxiliary.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/ChemicalAuxiliary.o ${OBJECTDIR}/src/core/chem/ChemicalAuxiliary_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/SimCoefCalculator_nomain.o: ${OBJECTDIR}/src/core/chem/SimCoefCalculator.o src/core/chem/SimCoefCalculator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/SimCoefCalculator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/SimCoefCalculator_nomain.o src/core/chem/SimCoefCalculator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/SimCoefCalculator.o ${OBJECTDIR}/src/core/chem/SimCoefCalculator_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr.o src/core/chem/fingerprintStrategy/AtomPairsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr_nomain.o src/core/chem/fingerprintStrategy/AtomPairsFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/AtomPairsFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy.o src/core/chem/fingerprintStrategy/FingerprintStrategy.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy_nomain.o src/core/chem/fingerprintStrategy/FingerprintStrategy.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/FingerprintStrategy_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr.o src/core/chem/fingerprintStrategy/MorganFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr_nomain.o src/core/chem/fingerprintStrategy/MorganFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/MorganFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.o src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1_nomain.o src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr1_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.o src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2_nomain.o src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolLayeredFngpr2_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr.o src/core/chem/fingerprintStrategy/TopolSingleFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr_nomain.o src/core/chem/fingerprintStrategy/TopolSingleFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolSingleFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr.o src/core/chem/fingerprintStrategy/TopolTorsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr_nomain.o src/core/chem/fingerprintStrategy/TopolTorsFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/TopolTorsFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr_nomain.o: ${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr.o src/core/chem/fingerprintStrategy/VectorFpFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr_nomain.o src/core/chem/fingerprintStrategy/VectorFpFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr.o ${OBJECTDIR}/src/core/chem/fingerprintStrategy/VectorFpFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphing/Morphing_nomain.o: ${OBJECTDIR}/src/core/chem/morphing/Morphing.o src/core/chem/morphing/Morphing.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphing/Morphing.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphing/Morphing_nomain.o src/core/chem/morphing/Morphing.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphing/Morphing.o ${OBJECTDIR}/src/core/chem/morphing/Morphing_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphing/MorphingData_nomain.o: ${OBJECTDIR}/src/core/chem/morphing/MorphingData.o src/core/chem/morphing/MorphingData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphing/MorphingData.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphing/MorphingData_nomain.o src/core/chem/morphing/MorphingData.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphing/MorphingData.o ${OBJECTDIR}/src/core/chem/morphing/MorphingData_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphing/MorphingFtors_nomain.o: ${OBJECTDIR}/src/core/chem/morphing/MorphingFtors.o src/core/chem/morphing/MorphingFtors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphing/MorphingFtors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphing/MorphingFtors_nomain.o src/core/chem/morphing/MorphingFtors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphing/MorphingFtors.o ${OBJECTDIR}/src/core/chem/morphing/MorphingFtors_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom.o src/core/chem/morphingStrategy/OpAddAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom_nomain.o src/core/chem/morphingStrategy/OpAddAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddAtom_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond.o src/core/chem/morphingStrategy/OpAddBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond_nomain.o src/core/chem/morphingStrategy/OpAddBond.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpAddBond_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction.o src/core/chem/morphingStrategy/OpBondContraction.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction_nomain.o src/core/chem/morphingStrategy/OpBondContraction.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondContraction_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute.o src/core/chem/morphingStrategy/OpBondReroute.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute_nomain.o src/core/chem/morphingStrategy/OpBondReroute.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpBondReroute_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom.o src/core/chem/morphingStrategy/OpInterlayAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom_nomain.o src/core/chem/morphingStrategy/OpInterlayAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpInterlayAtom_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom.o src/core/chem/morphingStrategy/OpMutateAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom_nomain.o src/core/chem/morphingStrategy/OpMutateAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpMutateAtom_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom.o src/core/chem/morphingStrategy/OpRemoveAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom_nomain.o src/core/chem/morphingStrategy/OpRemoveAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveAtom_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond_nomain.o: ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond.o src/core/chem/morphingStrategy/OpRemoveBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond_nomain.o src/core/chem/morphingStrategy/OpRemoveBond.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond.o ${OBJECTDIR}/src/core/chem/morphingStrategy/OpRemoveBond_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef.o src/core/chem/simCoefStrategy/AllBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef_nomain.o src/core/chem/simCoefStrategy/AllBitSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/AllBitSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef.o src/core/chem/simCoefStrategy/AsymmetricSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef_nomain.o src/core/chem/simCoefStrategy/AsymmetricSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/AsymmetricSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.o src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef_nomain.o src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/BraunBlanquetSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef.o src/core/chem/simCoefStrategy/CosineSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef_nomain.o src/core/chem/simCoefStrategy/CosineSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/CosineSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef.o src/core/chem/simCoefStrategy/DiceSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef_nomain.o src/core/chem/simCoefStrategy/DiceSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/DiceSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef.o src/core/chem/simCoefStrategy/KulczynskiSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef_nomain.o src/core/chem/simCoefStrategy/KulczynskiSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/KulczynskiSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef.o src/core/chem/simCoefStrategy/McConnaugheySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef_nomain.o src/core/chem/simCoefStrategy/McConnaugheySimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/McConnaugheySimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef.o src/core/chem/simCoefStrategy/OnBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef_nomain.o src/core/chem/simCoefStrategy/OnBitSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/OnBitSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef.o src/core/chem/simCoefStrategy/RusselSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef_nomain.o src/core/chem/simCoefStrategy/RusselSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/RusselSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef.o src/core/chem/simCoefStrategy/SokalSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef_nomain.o src/core/chem/simCoefStrategy/SokalSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/SokalSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef.o src/core/chem/simCoefStrategy/TanimotoSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef_nomain.o src/core/chem/simCoefStrategy/TanimotoSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/TanimotoSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef_nomain.o: ${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef.o src/core/chem/simCoefStrategy/TverskySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef_nomain.o src/core/chem/simCoefStrategy/TverskySimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef.o ${OBJECTDIR}/src/core/chem/simCoefStrategy/TverskySimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/SAScore_nomain.o: ${OBJECTDIR}/src/core/misc/SAScore.o src/core/misc/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/SAScore.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SAScore_nomain.o src/core/misc/SAScore.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/SAScore.o ${OBJECTDIR}/src/core/misc/SAScore_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/SAScore_data_loader_nomain.o: ${OBJECTDIR}/src/core/misc/SAScore_data_loader.o src/core/misc/SAScore_data_loader.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/SAScore_data_loader.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SAScore_data_loader_nomain.o src/core/misc/SAScore_data_loader.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/SAScore_data_loader.o ${OBJECTDIR}/src/core/misc/SAScore_data_loader_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/SynchRand_nomain.o: ${OBJECTDIR}/src/core/misc/SynchRand.o src/core/misc/SynchRand.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/SynchRand.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SynchRand_nomain.o src/core/misc/SynchRand.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/SynchRand.o ${OBJECTDIR}/src/core/misc/SynchRand_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/inout_nomain.o: ${OBJECTDIR}/src/core/misc/inout.o src/core/misc/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/inout.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/inout_nomain.o src/core/misc/inout.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/inout.o ${OBJECTDIR}/src/core/misc/inout_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/iteration_serializer_nomain.o: ${OBJECTDIR}/src/core/misc/iteration_serializer.o src/core/misc/iteration_serializer.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/iteration_serializer.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/iteration_serializer_nomain.o src/core/misc/iteration_serializer.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/iteration_serializer.o ${OBJECTDIR}/src/core/misc/iteration_serializer_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors_nomain.o: ${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors.o src/core/misc/selectors/chemoper_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors_nomain.o src/core/misc/selectors/chemoper_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors.o ${OBJECTDIR}/src/core/misc/selectors/chemoper_selectors_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors_nomain.o: ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o src/core/misc/selectors/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors_nomain.o src/core/misc/selectors/fingerprint_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors_nomain.o: ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o src/core/misc/selectors/simcoeff_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -DDBOOST_ALL_NO_LIB -DDBOOST_THREAD_USE_LIB -Isrc/ -Ideps/tbb/include/ -Iinclude/ -Ideps/rdkit/Code/ -Ideps/boost/ -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors_nomain.o src/core/misc/selectors/simcoeff_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} dist/lib/libmolpher.so

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
