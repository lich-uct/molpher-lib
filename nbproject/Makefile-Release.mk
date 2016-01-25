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
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/include/API.o \
	${OBJECTDIR}/include/operations/callbacks/callbacks.o \
	${OBJECTDIR}/include/operations/operations.o \
	${OBJECTDIR}/src/API/ExplorationParameters.o \
	${OBJECTDIR}/src/API/ExplorationTree.o \
	${OBJECTDIR}/src/API/ExplorationTreeSnapshot.o \
	${OBJECTDIR}/src/API/MolpherMol.o \
	${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback.o \
	${OBJECTDIR}/src/API/callbacks/TraverseCallback.o \
	${OBJECTDIR}/src/API/operations/ExtendTreeOper.o \
	${OBJECTDIR}/src/API/operations/FilterMorphsOper.o \
	${OBJECTDIR}/src/API/operations/FindLeavesOper.o \
	${OBJECTDIR}/src/API/operations/GenerateMorphsOper.o \
	${OBJECTDIR}/src/API/operations/PruneTreeOper.o \
	${OBJECTDIR}/src/API/operations/SortMorphsOper.o \
	${OBJECTDIR}/src/API/operations/TraverseOper.o \
	${OBJECTDIR}/src/API/operations/TreeOperation.o \
	${OBJECTDIR}/src/auxilliary/SAScore.o \
	${OBJECTDIR}/src/auxilliary/SynchRand.o \
	${OBJECTDIR}/src/auxilliary/inout.o \
	${OBJECTDIR}/src/auxilliary/iteration_serializer.o \
	${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors.o \
	${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors.o \
	${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors.o \
	${OBJECTDIR}/src/chem/ChemicalAuxiliary.o \
	${OBJECTDIR}/src/chem/SimCoefCalculator.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr.o \
	${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr.o \
	${OBJECTDIR}/src/chem/morphing/Morphing.o \
	${OBJECTDIR}/src/chem/morphing/MorphingData.o \
	${OBJECTDIR}/src/chem/morphing/MorphingFtors.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom.o \
	${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef.o \
	${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef.o \
	${OBJECTDIR}/src/swig/molpher.o \
	${OBJECTDIR}/src/swig/molpher_wrap.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1

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
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmolpher-lib.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/include/API.o: include/API.i 
	${MKDIR} -p ${OBJECTDIR}/include
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/include/API.o include/API.i

${OBJECTDIR}/include/operations/callbacks/callbacks.o: include/operations/callbacks/callbacks.i 
	${MKDIR} -p ${OBJECTDIR}/include/operations/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/include/operations/callbacks/callbacks.o include/operations/callbacks/callbacks.i

${OBJECTDIR}/include/operations/operations.o: include/operations/operations.i 
	${MKDIR} -p ${OBJECTDIR}/include/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/include/operations/operations.o include/operations/operations.i

${OBJECTDIR}/src/API/ExplorationParameters.o: src/API/ExplorationParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/ExplorationParameters.o src/API/ExplorationParameters.cpp

${OBJECTDIR}/src/API/ExplorationTree.o: src/API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/ExplorationTree.o src/API/ExplorationTree.cpp

${OBJECTDIR}/src/API/ExplorationTreeSnapshot.o: src/API/ExplorationTreeSnapshot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/ExplorationTreeSnapshot.o src/API/ExplorationTreeSnapshot.cpp

${OBJECTDIR}/src/API/MolpherMol.o: src/API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/MolpherMol.o src/API/MolpherMol.cpp

${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback.o: src/API/callbacks/EraseSubtreeCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback.o src/API/callbacks/EraseSubtreeCallback.cpp

${OBJECTDIR}/src/API/callbacks/TraverseCallback.o: src/API/callbacks/TraverseCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/callbacks/TraverseCallback.o src/API/callbacks/TraverseCallback.cpp

${OBJECTDIR}/src/API/operations/ExtendTreeOper.o: src/API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/ExtendTreeOper.o src/API/operations/ExtendTreeOper.cpp

${OBJECTDIR}/src/API/operations/FilterMorphsOper.o: src/API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/FilterMorphsOper.o src/API/operations/FilterMorphsOper.cpp

${OBJECTDIR}/src/API/operations/FindLeavesOper.o: src/API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/FindLeavesOper.o src/API/operations/FindLeavesOper.cpp

${OBJECTDIR}/src/API/operations/GenerateMorphsOper.o: src/API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/GenerateMorphsOper.o src/API/operations/GenerateMorphsOper.cpp

${OBJECTDIR}/src/API/operations/PruneTreeOper.o: src/API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/PruneTreeOper.o src/API/operations/PruneTreeOper.cpp

${OBJECTDIR}/src/API/operations/SortMorphsOper.o: src/API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/SortMorphsOper.o src/API/operations/SortMorphsOper.cpp

${OBJECTDIR}/src/API/operations/TraverseOper.o: src/API/operations/TraverseOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/TraverseOper.o src/API/operations/TraverseOper.cpp

${OBJECTDIR}/src/API/operations/TreeOperation.o: src/API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/TreeOperation.o src/API/operations/TreeOperation.cpp

${OBJECTDIR}/src/auxilliary/SAScore.o: src/auxilliary/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/SAScore.o src/auxilliary/SAScore.cpp

${OBJECTDIR}/src/auxilliary/SynchRand.o: src/auxilliary/SynchRand.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/SynchRand.o src/auxilliary/SynchRand.cpp

${OBJECTDIR}/src/auxilliary/inout.o: src/auxilliary/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/inout.o src/auxilliary/inout.cpp

${OBJECTDIR}/src/auxilliary/iteration_serializer.o: src/auxilliary/iteration_serializer.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/iteration_serializer.o src/auxilliary/iteration_serializer.cpp

${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors.o: src/auxilliary/selectors/chemoper_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors.o src/auxilliary/selectors/chemoper_selectors.cpp

${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors.o: src/auxilliary/selectors/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors.o src/auxilliary/selectors/fingerprint_selectors.cpp

${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors.o: src/auxilliary/selectors/simcoeff_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors.o src/auxilliary/selectors/simcoeff_selectors.cpp

${OBJECTDIR}/src/chem/ChemicalAuxiliary.o: src/chem/ChemicalAuxiliary.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/ChemicalAuxiliary.o src/chem/ChemicalAuxiliary.cpp

${OBJECTDIR}/src/chem/SimCoefCalculator.o: src/chem/SimCoefCalculator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/SimCoefCalculator.o src/chem/SimCoefCalculator.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr.o: src/chem/fingerprintStrategy/AtomPairsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr.o src/chem/fingerprintStrategy/AtomPairsFngpr.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy.o: src/chem/fingerprintStrategy/FingerprintStrategy.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy.o src/chem/fingerprintStrategy/FingerprintStrategy.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr.o: src/chem/fingerprintStrategy/MorganFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr.o src/chem/fingerprintStrategy/MorganFngpr.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1.o: src/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1.o src/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2.o: src/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2.o src/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr.o: src/chem/fingerprintStrategy/TopolSingleFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr.o src/chem/fingerprintStrategy/TopolSingleFngpr.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr.o: src/chem/fingerprintStrategy/TopolTorsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr.o src/chem/fingerprintStrategy/TopolTorsFngpr.cpp

${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr.o: src/chem/fingerprintStrategy/VectorFpFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr.o src/chem/fingerprintStrategy/VectorFpFngpr.cpp

${OBJECTDIR}/src/chem/morphing/Morphing.o: src/chem/morphing/Morphing.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphing/Morphing.o src/chem/morphing/Morphing.cpp

${OBJECTDIR}/src/chem/morphing/MorphingData.o: src/chem/morphing/MorphingData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphing/MorphingData.o src/chem/morphing/MorphingData.cpp

${OBJECTDIR}/src/chem/morphing/MorphingFtors.o: src/chem/morphing/MorphingFtors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphing/MorphingFtors.o src/chem/morphing/MorphingFtors.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom.o: src/chem/morphingStrategy/OpAddAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom.o src/chem/morphingStrategy/OpAddAtom.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond.o: src/chem/morphingStrategy/OpAddBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond.o src/chem/morphingStrategy/OpAddBond.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction.o: src/chem/morphingStrategy/OpBondContraction.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction.o src/chem/morphingStrategy/OpBondContraction.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute.o: src/chem/morphingStrategy/OpBondReroute.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute.o src/chem/morphingStrategy/OpBondReroute.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom.o: src/chem/morphingStrategy/OpInterlayAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom.o src/chem/morphingStrategy/OpInterlayAtom.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom.o: src/chem/morphingStrategy/OpMutateAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom.o src/chem/morphingStrategy/OpMutateAtom.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom.o: src/chem/morphingStrategy/OpRemoveAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom.o src/chem/morphingStrategy/OpRemoveAtom.cpp

${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond.o: src/chem/morphingStrategy/OpRemoveBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond.o src/chem/morphingStrategy/OpRemoveBond.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef.o: src/chem/simCoefStrategy/AllBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef.o src/chem/simCoefStrategy/AllBitSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef.o: src/chem/simCoefStrategy/AsymmetricSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef.o src/chem/simCoefStrategy/AsymmetricSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef.o: src/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef.o src/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef.o: src/chem/simCoefStrategy/CosineSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef.o src/chem/simCoefStrategy/CosineSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef.o: src/chem/simCoefStrategy/DiceSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef.o src/chem/simCoefStrategy/DiceSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef.o: src/chem/simCoefStrategy/KulczynskiSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef.o src/chem/simCoefStrategy/KulczynskiSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef.o: src/chem/simCoefStrategy/McConnaugheySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef.o src/chem/simCoefStrategy/McConnaugheySimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef.o: src/chem/simCoefStrategy/OnBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef.o src/chem/simCoefStrategy/OnBitSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef.o: src/chem/simCoefStrategy/RusselSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef.o src/chem/simCoefStrategy/RusselSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef.o: src/chem/simCoefStrategy/SokalSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef.o src/chem/simCoefStrategy/SokalSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef.o: src/chem/simCoefStrategy/TanimotoSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef.o src/chem/simCoefStrategy/TanimotoSimCoef.cpp

${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef.o: src/chem/simCoefStrategy/TverskySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef.o src/chem/simCoefStrategy/TverskySimCoef.cpp

${OBJECTDIR}/src/swig/molpher.o: src/swig/molpher.i 
	${MKDIR} -p ${OBJECTDIR}/src/swig
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/swig/molpher.o src/swig/molpher.i

${OBJECTDIR}/src/swig/molpher_wrap.o: src/swig/molpher_wrap.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/swig
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/swig/molpher_wrap.o src/swig/molpher_wrap.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/APITests.o ${TESTDIR}/tests/testrunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} `cppunit-config --libs`   


${TESTDIR}/tests/APITests.o: tests/APITests.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/APITests.o tests/APITests.cpp


${TESTDIR}/tests/testrunner.o: tests/testrunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/testrunner.o tests/testrunner.cpp


${OBJECTDIR}/include/API_nomain.o: ${OBJECTDIR}/include/API.o include/API.i 
	${MKDIR} -p ${OBJECTDIR}/include
	@NMOUTPUT=`${NM} ${OBJECTDIR}/include/API.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/include/API_nomain.o include/API.i;\
	else  \
	    ${CP} ${OBJECTDIR}/include/API.o ${OBJECTDIR}/include/API_nomain.o;\
	fi

${OBJECTDIR}/include/operations/callbacks/callbacks_nomain.o: ${OBJECTDIR}/include/operations/callbacks/callbacks.o include/operations/callbacks/callbacks.i 
	${MKDIR} -p ${OBJECTDIR}/include/operations/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/include/operations/callbacks/callbacks.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/include/operations/callbacks/callbacks_nomain.o include/operations/callbacks/callbacks.i;\
	else  \
	    ${CP} ${OBJECTDIR}/include/operations/callbacks/callbacks.o ${OBJECTDIR}/include/operations/callbacks/callbacks_nomain.o;\
	fi

${OBJECTDIR}/include/operations/operations_nomain.o: ${OBJECTDIR}/include/operations/operations.o include/operations/operations.i 
	${MKDIR} -p ${OBJECTDIR}/include/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/include/operations/operations.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/include/operations/operations_nomain.o include/operations/operations.i;\
	else  \
	    ${CP} ${OBJECTDIR}/include/operations/operations.o ${OBJECTDIR}/include/operations/operations_nomain.o;\
	fi

${OBJECTDIR}/src/API/ExplorationParameters_nomain.o: ${OBJECTDIR}/src/API/ExplorationParameters.o src/API/ExplorationParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/ExplorationParameters.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/ExplorationParameters_nomain.o src/API/ExplorationParameters.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/ExplorationParameters.o ${OBJECTDIR}/src/API/ExplorationParameters_nomain.o;\
	fi

${OBJECTDIR}/src/API/ExplorationTree_nomain.o: ${OBJECTDIR}/src/API/ExplorationTree.o src/API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/ExplorationTree.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/ExplorationTree_nomain.o src/API/ExplorationTree.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/ExplorationTree.o ${OBJECTDIR}/src/API/ExplorationTree_nomain.o;\
	fi

${OBJECTDIR}/src/API/ExplorationTreeSnapshot_nomain.o: ${OBJECTDIR}/src/API/ExplorationTreeSnapshot.o src/API/ExplorationTreeSnapshot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/ExplorationTreeSnapshot.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/ExplorationTreeSnapshot_nomain.o src/API/ExplorationTreeSnapshot.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/ExplorationTreeSnapshot.o ${OBJECTDIR}/src/API/ExplorationTreeSnapshot_nomain.o;\
	fi

${OBJECTDIR}/src/API/MolpherMol_nomain.o: ${OBJECTDIR}/src/API/MolpherMol.o src/API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/MolpherMol.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/MolpherMol_nomain.o src/API/MolpherMol.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/MolpherMol.o ${OBJECTDIR}/src/API/MolpherMol_nomain.o;\
	fi

${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback_nomain.o: ${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback.o src/API/callbacks/EraseSubtreeCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback_nomain.o src/API/callbacks/EraseSubtreeCallback.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback.o ${OBJECTDIR}/src/API/callbacks/EraseSubtreeCallback_nomain.o;\
	fi

${OBJECTDIR}/src/API/callbacks/TraverseCallback_nomain.o: ${OBJECTDIR}/src/API/callbacks/TraverseCallback.o src/API/callbacks/TraverseCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/callbacks/TraverseCallback.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/callbacks/TraverseCallback_nomain.o src/API/callbacks/TraverseCallback.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/callbacks/TraverseCallback.o ${OBJECTDIR}/src/API/callbacks/TraverseCallback_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/ExtendTreeOper_nomain.o: ${OBJECTDIR}/src/API/operations/ExtendTreeOper.o src/API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/ExtendTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/ExtendTreeOper_nomain.o src/API/operations/ExtendTreeOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/ExtendTreeOper.o ${OBJECTDIR}/src/API/operations/ExtendTreeOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/FilterMorphsOper_nomain.o: ${OBJECTDIR}/src/API/operations/FilterMorphsOper.o src/API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/FilterMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/FilterMorphsOper_nomain.o src/API/operations/FilterMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/FilterMorphsOper.o ${OBJECTDIR}/src/API/operations/FilterMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/FindLeavesOper_nomain.o: ${OBJECTDIR}/src/API/operations/FindLeavesOper.o src/API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/FindLeavesOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/FindLeavesOper_nomain.o src/API/operations/FindLeavesOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/FindLeavesOper.o ${OBJECTDIR}/src/API/operations/FindLeavesOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/GenerateMorphsOper_nomain.o: ${OBJECTDIR}/src/API/operations/GenerateMorphsOper.o src/API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/GenerateMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/GenerateMorphsOper_nomain.o src/API/operations/GenerateMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/GenerateMorphsOper.o ${OBJECTDIR}/src/API/operations/GenerateMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/PruneTreeOper_nomain.o: ${OBJECTDIR}/src/API/operations/PruneTreeOper.o src/API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/PruneTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/PruneTreeOper_nomain.o src/API/operations/PruneTreeOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/PruneTreeOper.o ${OBJECTDIR}/src/API/operations/PruneTreeOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/SortMorphsOper_nomain.o: ${OBJECTDIR}/src/API/operations/SortMorphsOper.o src/API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/SortMorphsOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/SortMorphsOper_nomain.o src/API/operations/SortMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/SortMorphsOper.o ${OBJECTDIR}/src/API/operations/SortMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/TraverseOper_nomain.o: ${OBJECTDIR}/src/API/operations/TraverseOper.o src/API/operations/TraverseOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/TraverseOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/TraverseOper_nomain.o src/API/operations/TraverseOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/TraverseOper.o ${OBJECTDIR}/src/API/operations/TraverseOper_nomain.o;\
	fi

${OBJECTDIR}/src/API/operations/TreeOperation_nomain.o: ${OBJECTDIR}/src/API/operations/TreeOperation.o src/API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/API/operations/TreeOperation.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/API/operations/TreeOperation_nomain.o src/API/operations/TreeOperation.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/API/operations/TreeOperation.o ${OBJECTDIR}/src/API/operations/TreeOperation_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/SAScore_nomain.o: ${OBJECTDIR}/src/auxilliary/SAScore.o src/auxilliary/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/SAScore.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/SAScore_nomain.o src/auxilliary/SAScore.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/SAScore.o ${OBJECTDIR}/src/auxilliary/SAScore_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/SynchRand_nomain.o: ${OBJECTDIR}/src/auxilliary/SynchRand.o src/auxilliary/SynchRand.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/SynchRand.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/SynchRand_nomain.o src/auxilliary/SynchRand.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/SynchRand.o ${OBJECTDIR}/src/auxilliary/SynchRand_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/inout_nomain.o: ${OBJECTDIR}/src/auxilliary/inout.o src/auxilliary/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/inout.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/inout_nomain.o src/auxilliary/inout.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/inout.o ${OBJECTDIR}/src/auxilliary/inout_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/iteration_serializer_nomain.o: ${OBJECTDIR}/src/auxilliary/iteration_serializer.o src/auxilliary/iteration_serializer.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/iteration_serializer.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/iteration_serializer_nomain.o src/auxilliary/iteration_serializer.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/iteration_serializer.o ${OBJECTDIR}/src/auxilliary/iteration_serializer_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors_nomain.o: ${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors.o src/auxilliary/selectors/chemoper_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors_nomain.o src/auxilliary/selectors/chemoper_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors.o ${OBJECTDIR}/src/auxilliary/selectors/chemoper_selectors_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors_nomain.o: ${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors.o src/auxilliary/selectors/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors_nomain.o src/auxilliary/selectors/fingerprint_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors.o ${OBJECTDIR}/src/auxilliary/selectors/fingerprint_selectors_nomain.o;\
	fi

${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors_nomain.o: ${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors.o src/auxilliary/selectors/simcoeff_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/auxilliary/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors_nomain.o src/auxilliary/selectors/simcoeff_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors.o ${OBJECTDIR}/src/auxilliary/selectors/simcoeff_selectors_nomain.o;\
	fi

${OBJECTDIR}/src/chem/ChemicalAuxiliary_nomain.o: ${OBJECTDIR}/src/chem/ChemicalAuxiliary.o src/chem/ChemicalAuxiliary.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/ChemicalAuxiliary.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/ChemicalAuxiliary_nomain.o src/chem/ChemicalAuxiliary.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/ChemicalAuxiliary.o ${OBJECTDIR}/src/chem/ChemicalAuxiliary_nomain.o;\
	fi

${OBJECTDIR}/src/chem/SimCoefCalculator_nomain.o: ${OBJECTDIR}/src/chem/SimCoefCalculator.o src/chem/SimCoefCalculator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/SimCoefCalculator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/SimCoefCalculator_nomain.o src/chem/SimCoefCalculator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/SimCoefCalculator.o ${OBJECTDIR}/src/chem/SimCoefCalculator_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr.o src/chem/fingerprintStrategy/AtomPairsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr_nomain.o src/chem/fingerprintStrategy/AtomPairsFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr.o ${OBJECTDIR}/src/chem/fingerprintStrategy/AtomPairsFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy.o src/chem/fingerprintStrategy/FingerprintStrategy.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy_nomain.o src/chem/fingerprintStrategy/FingerprintStrategy.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy.o ${OBJECTDIR}/src/chem/fingerprintStrategy/FingerprintStrategy_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr.o src/chem/fingerprintStrategy/MorganFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr_nomain.o src/chem/fingerprintStrategy/MorganFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr.o ${OBJECTDIR}/src/chem/fingerprintStrategy/MorganFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1.o src/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1_nomain.o src/chem/fingerprintStrategy/TopolLayeredFngpr1.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1.o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr1_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2.o src/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2_nomain.o src/chem/fingerprintStrategy/TopolLayeredFngpr2.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2.o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolLayeredFngpr2_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr.o src/chem/fingerprintStrategy/TopolSingleFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr_nomain.o src/chem/fingerprintStrategy/TopolSingleFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr.o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolSingleFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr.o src/chem/fingerprintStrategy/TopolTorsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr_nomain.o src/chem/fingerprintStrategy/TopolTorsFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr.o ${OBJECTDIR}/src/chem/fingerprintStrategy/TopolTorsFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr_nomain.o: ${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr.o src/chem/fingerprintStrategy/VectorFpFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/fingerprintStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr_nomain.o src/chem/fingerprintStrategy/VectorFpFngpr.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr.o ${OBJECTDIR}/src/chem/fingerprintStrategy/VectorFpFngpr_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphing/Morphing_nomain.o: ${OBJECTDIR}/src/chem/morphing/Morphing.o src/chem/morphing/Morphing.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphing/Morphing.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphing/Morphing_nomain.o src/chem/morphing/Morphing.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphing/Morphing.o ${OBJECTDIR}/src/chem/morphing/Morphing_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphing/MorphingData_nomain.o: ${OBJECTDIR}/src/chem/morphing/MorphingData.o src/chem/morphing/MorphingData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphing/MorphingData.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphing/MorphingData_nomain.o src/chem/morphing/MorphingData.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphing/MorphingData.o ${OBJECTDIR}/src/chem/morphing/MorphingData_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphing/MorphingFtors_nomain.o: ${OBJECTDIR}/src/chem/morphing/MorphingFtors.o src/chem/morphing/MorphingFtors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphing/MorphingFtors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphing/MorphingFtors_nomain.o src/chem/morphing/MorphingFtors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphing/MorphingFtors.o ${OBJECTDIR}/src/chem/morphing/MorphingFtors_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom.o src/chem/morphingStrategy/OpAddAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom_nomain.o src/chem/morphingStrategy/OpAddAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom.o ${OBJECTDIR}/src/chem/morphingStrategy/OpAddAtom_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond.o src/chem/morphingStrategy/OpAddBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond_nomain.o src/chem/morphingStrategy/OpAddBond.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond.o ${OBJECTDIR}/src/chem/morphingStrategy/OpAddBond_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction.o src/chem/morphingStrategy/OpBondContraction.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction_nomain.o src/chem/morphingStrategy/OpBondContraction.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction.o ${OBJECTDIR}/src/chem/morphingStrategy/OpBondContraction_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute.o src/chem/morphingStrategy/OpBondReroute.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute_nomain.o src/chem/morphingStrategy/OpBondReroute.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute.o ${OBJECTDIR}/src/chem/morphingStrategy/OpBondReroute_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom.o src/chem/morphingStrategy/OpInterlayAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom_nomain.o src/chem/morphingStrategy/OpInterlayAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom.o ${OBJECTDIR}/src/chem/morphingStrategy/OpInterlayAtom_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom.o src/chem/morphingStrategy/OpMutateAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom_nomain.o src/chem/morphingStrategy/OpMutateAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom.o ${OBJECTDIR}/src/chem/morphingStrategy/OpMutateAtom_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom.o src/chem/morphingStrategy/OpRemoveAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom_nomain.o src/chem/morphingStrategy/OpRemoveAtom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom.o ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveAtom_nomain.o;\
	fi

${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond_nomain.o: ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond.o src/chem/morphingStrategy/OpRemoveBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/morphingStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond_nomain.o src/chem/morphingStrategy/OpRemoveBond.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond.o ${OBJECTDIR}/src/chem/morphingStrategy/OpRemoveBond_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef.o src/chem/simCoefStrategy/AllBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef_nomain.o src/chem/simCoefStrategy/AllBitSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/AllBitSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef.o src/chem/simCoefStrategy/AsymmetricSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef_nomain.o src/chem/simCoefStrategy/AsymmetricSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/AsymmetricSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef.o src/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef_nomain.o src/chem/simCoefStrategy/BraunBlanquetSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/BraunBlanquetSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef.o src/chem/simCoefStrategy/CosineSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef_nomain.o src/chem/simCoefStrategy/CosineSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/CosineSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef.o src/chem/simCoefStrategy/DiceSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef_nomain.o src/chem/simCoefStrategy/DiceSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/DiceSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef.o src/chem/simCoefStrategy/KulczynskiSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef_nomain.o src/chem/simCoefStrategy/KulczynskiSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/KulczynskiSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef.o src/chem/simCoefStrategy/McConnaugheySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef_nomain.o src/chem/simCoefStrategy/McConnaugheySimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/McConnaugheySimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef.o src/chem/simCoefStrategy/OnBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef_nomain.o src/chem/simCoefStrategy/OnBitSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/OnBitSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef.o src/chem/simCoefStrategy/RusselSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef_nomain.o src/chem/simCoefStrategy/RusselSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/RusselSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef.o src/chem/simCoefStrategy/SokalSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef_nomain.o src/chem/simCoefStrategy/SokalSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/SokalSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef.o src/chem/simCoefStrategy/TanimotoSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef_nomain.o src/chem/simCoefStrategy/TanimotoSimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/TanimotoSimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef_nomain.o: ${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef.o src/chem/simCoefStrategy/TverskySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/chem/simCoefStrategy
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef_nomain.o src/chem/simCoefStrategy/TverskySimCoef.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef.o ${OBJECTDIR}/src/chem/simCoefStrategy/TverskySimCoef_nomain.o;\
	fi

${OBJECTDIR}/src/swig/molpher_nomain.o: ${OBJECTDIR}/src/swig/molpher.o src/swig/molpher.i 
	${MKDIR} -p ${OBJECTDIR}/src/swig
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/swig/molpher.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/swig/molpher_nomain.o src/swig/molpher.i;\
	else  \
	    ${CP} ${OBJECTDIR}/src/swig/molpher.o ${OBJECTDIR}/src/swig/molpher_nomain.o;\
	fi

${OBJECTDIR}/src/swig/molpher_wrap_nomain.o: ${OBJECTDIR}/src/swig/molpher_wrap.o src/swig/molpher_wrap.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/swig
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/swig/molpher_wrap.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/swig/molpher_wrap_nomain.o src/swig/molpher_wrap.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/swig/molpher_wrap.o ${OBJECTDIR}/src/swig/molpher_wrap_nomain.o;\
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
