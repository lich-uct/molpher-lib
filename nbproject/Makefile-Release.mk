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
	${OBJECTDIR}/src/core/API/ExplorationParameters.o \
	${OBJECTDIR}/src/core/API/ExplorationTree.o \
	${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot.o \
	${OBJECTDIR}/src/core/API/MolpherMol.o \
	${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback.o \
	${OBJECTDIR}/src/core/API/callbacks/TraverseCallback.o \
	${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o \
	${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o \
	${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o \
	${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o \
	${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o \
	${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o \
	${OBJECTDIR}/src/core/API/operations/TraverseOper.o \
	${OBJECTDIR}/src/core/API/operations/TreeOperation.o \
	${OBJECTDIR}/src/core/misc/SAScore.o \
	${OBJECTDIR}/src/core/misc/inout.o \
	${OBJECTDIR}/src/core/misc/iteration_serializer.o \
	${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o \
	${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o \
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

${OBJECTDIR}/src/core/API/ExplorationParameters.o: src/core/API/ExplorationParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationParameters.o src/core/API/ExplorationParameters.cpp

${OBJECTDIR}/src/core/API/ExplorationTree.o: src/core/API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationTree.o src/core/API/ExplorationTree.cpp

${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot.o: src/core/API/ExplorationTreeSnapshot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot.o src/core/API/ExplorationTreeSnapshot.cpp

${OBJECTDIR}/src/core/API/MolpherMol.o: src/core/API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/MolpherMol.o src/core/API/MolpherMol.cpp

${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback.o: src/core/API/callbacks/EraseSubtreeCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback.o src/core/API/callbacks/EraseSubtreeCallback.cpp

${OBJECTDIR}/src/core/API/callbacks/TraverseCallback.o: src/core/API/callbacks/TraverseCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/callbacks
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/callbacks/TraverseCallback.o src/core/API/callbacks/TraverseCallback.cpp

${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o: src/core/API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o src/core/API/operations/ExtendTreeOper.cpp

${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o: src/core/API/operations/FilterMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper.o src/core/API/operations/FilterMorphsOper.cpp

${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o: src/core/API/operations/FindLeavesOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FindLeavesOper.o src/core/API/operations/FindLeavesOper.cpp

${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o: src/core/API/operations/GenerateMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper.o src/core/API/operations/GenerateMorphsOper.cpp

${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o: src/core/API/operations/PruneTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/PruneTreeOper.o src/core/API/operations/PruneTreeOper.cpp

${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o: src/core/API/operations/SortMorphsOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o src/core/API/operations/SortMorphsOper.cpp

${OBJECTDIR}/src/core/API/operations/TraverseOper.o: src/core/API/operations/TraverseOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/TraverseOper.o src/core/API/operations/TraverseOper.cpp

${OBJECTDIR}/src/core/API/operations/TreeOperation.o: src/core/API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/TreeOperation.o src/core/API/operations/TreeOperation.cpp

${OBJECTDIR}/src/core/misc/SAScore.o: src/core/misc/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SAScore.o src/core/misc/SAScore.cpp

${OBJECTDIR}/src/core/misc/inout.o: src/core/misc/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/inout.o src/core/misc/inout.cpp

${OBJECTDIR}/src/core/misc/iteration_serializer.o: src/core/misc/iteration_serializer.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/iteration_serializer.o src/core/misc/iteration_serializer.cpp

${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o: src/core/misc/selectors/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o src/core/misc/selectors/fingerprint_selectors.cpp

${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o: src/core/misc/selectors/simcoeff_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o src/core/misc/selectors/simcoeff_selectors.cpp

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

${OBJECTDIR}/src/core/API/ExplorationParameters_nomain.o: ${OBJECTDIR}/src/core/API/ExplorationParameters.o src/core/API/ExplorationParameters.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/ExplorationParameters.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationParameters_nomain.o src/core/API/ExplorationParameters.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/ExplorationParameters.o ${OBJECTDIR}/src/core/API/ExplorationParameters_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/ExplorationTree_nomain.o: ${OBJECTDIR}/src/core/API/ExplorationTree.o src/core/API/ExplorationTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/ExplorationTree.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationTree_nomain.o src/core/API/ExplorationTree.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/ExplorationTree.o ${OBJECTDIR}/src/core/API/ExplorationTree_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot_nomain.o: ${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot.o src/core/API/ExplorationTreeSnapshot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot_nomain.o src/core/API/ExplorationTreeSnapshot.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot.o ${OBJECTDIR}/src/core/API/ExplorationTreeSnapshot_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/MolpherMol_nomain.o: ${OBJECTDIR}/src/core/API/MolpherMol.o src/core/API/MolpherMol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/MolpherMol.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/MolpherMol_nomain.o src/core/API/MolpherMol.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/MolpherMol.o ${OBJECTDIR}/src/core/API/MolpherMol_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback_nomain.o: ${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback.o src/core/API/callbacks/EraseSubtreeCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback_nomain.o src/core/API/callbacks/EraseSubtreeCallback.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback.o ${OBJECTDIR}/src/core/API/callbacks/EraseSubtreeCallback_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/callbacks/TraverseCallback_nomain.o: ${OBJECTDIR}/src/core/API/callbacks/TraverseCallback.o src/core/API/callbacks/TraverseCallback.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/callbacks
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/callbacks/TraverseCallback.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/callbacks/TraverseCallback_nomain.o src/core/API/callbacks/TraverseCallback.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/callbacks/TraverseCallback.o ${OBJECTDIR}/src/core/API/callbacks/TraverseCallback_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/ExtendTreeOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o src/core/API/operations/ExtendTreeOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/ExtendTreeOper_nomain.o src/core/API/operations/ExtendTreeOper.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FilterMorphsOper_nomain.o src/core/API/operations/FilterMorphsOper.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/FindLeavesOper_nomain.o src/core/API/operations/FindLeavesOper.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/GenerateMorphsOper_nomain.o src/core/API/operations/GenerateMorphsOper.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/PruneTreeOper_nomain.o src/core/API/operations/PruneTreeOper.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/SortMorphsOper_nomain.o src/core/API/operations/SortMorphsOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/SortMorphsOper.o ${OBJECTDIR}/src/core/API/operations/SortMorphsOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/TraverseOper_nomain.o: ${OBJECTDIR}/src/core/API/operations/TraverseOper.o src/core/API/operations/TraverseOper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/TraverseOper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/TraverseOper_nomain.o src/core/API/operations/TraverseOper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/TraverseOper.o ${OBJECTDIR}/src/core/API/operations/TraverseOper_nomain.o;\
	fi

${OBJECTDIR}/src/core/API/operations/TreeOperation_nomain.o: ${OBJECTDIR}/src/core/API/operations/TreeOperation.o src/core/API/operations/TreeOperation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/API/operations
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/API/operations/TreeOperation.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/API/operations/TreeOperation_nomain.o src/core/API/operations/TreeOperation.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/API/operations/TreeOperation.o ${OBJECTDIR}/src/core/API/operations/TreeOperation_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/SAScore_nomain.o: ${OBJECTDIR}/src/core/misc/SAScore.o src/core/misc/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/SAScore.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/SAScore_nomain.o src/core/misc/SAScore.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/SAScore.o ${OBJECTDIR}/src/core/misc/SAScore_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/inout_nomain.o: ${OBJECTDIR}/src/core/misc/inout.o src/core/misc/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/inout.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/inout_nomain.o src/core/misc/inout.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/iteration_serializer_nomain.o src/core/misc/iteration_serializer.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/iteration_serializer.o ${OBJECTDIR}/src/core/misc/iteration_serializer_nomain.o;\
	fi

${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors_nomain.o: ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o src/core/misc/selectors/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/core/misc/selectors
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/fingerprint_selectors_nomain.o src/core/misc/selectors/fingerprint_selectors.cpp;\
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
	    $(COMPILE.cc) -O2 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors_nomain.o src/core/misc/selectors/simcoeff_selectors.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors.o ${OBJECTDIR}/src/core/misc/selectors/simcoeff_selectors_nomain.o;\
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
