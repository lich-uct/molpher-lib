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
CC=gcc-4.6
CCC=g++-4.6
CXX=g++-4.6
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU_4_6-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Linux64_Debug_library
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1270477542/chemoper_selectors.o \
	${OBJECTDIR}/_ext/1270477542/dimred_selectors.o \
	${OBJECTDIR}/_ext/1270477542/fingerprint_selectors.o \
	${OBJECTDIR}/_ext/1270477542/inout.o \
	${OBJECTDIR}/_ext/1270477542/iteration_serializer.o \
	${OBJECTDIR}/_ext/1270477542/simcoeff_selectors.o \
	${OBJECTDIR}/BackendCommunicator.o \
	${OBJECTDIR}/auxiliary/SynchRand.o \
	${OBJECTDIR}/chem/ChemicalAuxiliary.o \
	${OBJECTDIR}/chem/SimCoefCalculator.o \
	${OBJECTDIR}/chem/fingerprintStrategy/AtomPairsFngpr.o \
	${OBJECTDIR}/chem/fingerprintStrategy/FingerprintStrategy.o \
	${OBJECTDIR}/chem/fingerprintStrategy/MorganFngpr.o \
	${OBJECTDIR}/chem/fingerprintStrategy/TopolLayeredFngpr1.o \
	${OBJECTDIR}/chem/fingerprintStrategy/TopolLayeredFngpr2.o \
	${OBJECTDIR}/chem/fingerprintStrategy/TopolSingleFngpr.o \
	${OBJECTDIR}/chem/fingerprintStrategy/TopolTorsFngpr.o \
	${OBJECTDIR}/chem/fingerprintStrategy/VectorFpFngpr.o \
	${OBJECTDIR}/chem/morphing/Morphing.o \
	${OBJECTDIR}/chem/morphing/MorphingData.o \
	${OBJECTDIR}/chem/morphing/MorphingFtors.o \
	${OBJECTDIR}/chem/morphingStrategy/OpAddAtom.o \
	${OBJECTDIR}/chem/morphingStrategy/OpAddBond.o \
	${OBJECTDIR}/chem/morphingStrategy/OpBondContraction.o \
	${OBJECTDIR}/chem/morphingStrategy/OpBondReroute.o \
	${OBJECTDIR}/chem/morphingStrategy/OpInterlayAtom.o \
	${OBJECTDIR}/chem/morphingStrategy/OpMutateAtom.o \
	${OBJECTDIR}/chem/morphingStrategy/OpRemoveAtom.o \
	${OBJECTDIR}/chem/morphingStrategy/OpRemoveBond.o \
	${OBJECTDIR}/chem/simCoefStrategy/AllBitSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/AsymmetricSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/BraunBlanquetSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/CosineSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/DiceSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/KulczynskiSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/McConnaugheySimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/OnBitSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/RusselSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/SokalSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/TanimotoSimCoef.o \
	${OBJECTDIR}/chem/simCoefStrategy/TverskySimCoef.o \
	${OBJECTDIR}/coord/KamadaKawaiReducer.o \
	${OBJECTDIR}/coord/PcaReducer.o \
	${OBJECTDIR}/coord/ReducerFactory.o \
	${OBJECTDIR}/core/JobManager.o \
	${OBJECTDIR}/core/NeighborhoodGenerator.o \
	${OBJECTDIR}/core/NeighborhoodTaskQueue.o \
	${OBJECTDIR}/core/PathFinder.o \
	${OBJECTDIR}/core/PathFinderContext.o \
	${OBJECTDIR}/extensions/SAScore.o \
	${OBJECTDIR}/tests/MorphingTest.o


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
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libbackend.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libbackend.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libbackend.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libbackend.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libbackend.a

${OBJECTDIR}/_ext/1270477542/chemoper_selectors.o: ../common/chemoper_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1270477542
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1270477542/chemoper_selectors.o ../common/chemoper_selectors.cpp

${OBJECTDIR}/_ext/1270477542/dimred_selectors.o: ../common/dimred_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1270477542
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1270477542/dimred_selectors.o ../common/dimred_selectors.cpp

${OBJECTDIR}/_ext/1270477542/fingerprint_selectors.o: ../common/fingerprint_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1270477542
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1270477542/fingerprint_selectors.o ../common/fingerprint_selectors.cpp

${OBJECTDIR}/_ext/1270477542/inout.o: ../common/inout.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1270477542
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1270477542/inout.o ../common/inout.cpp

${OBJECTDIR}/_ext/1270477542/iteration_serializer.o: ../common/iteration_serializer.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1270477542
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1270477542/iteration_serializer.o ../common/iteration_serializer.cpp

${OBJECTDIR}/_ext/1270477542/simcoeff_selectors.o: ../common/simcoeff_selectors.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1270477542
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1270477542/simcoeff_selectors.o ../common/simcoeff_selectors.cpp

${OBJECTDIR}/BackendCommunicator.o: BackendCommunicator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BackendCommunicator.o BackendCommunicator.cpp

${OBJECTDIR}/auxiliary/SynchRand.o: auxiliary/SynchRand.cpp 
	${MKDIR} -p ${OBJECTDIR}/auxiliary
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/auxiliary/SynchRand.o auxiliary/SynchRand.cpp

${OBJECTDIR}/chem/ChemicalAuxiliary.o: chem/ChemicalAuxiliary.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/ChemicalAuxiliary.o chem/ChemicalAuxiliary.cpp

${OBJECTDIR}/chem/SimCoefCalculator.o: chem/SimCoefCalculator.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/SimCoefCalculator.o chem/SimCoefCalculator.cpp

${OBJECTDIR}/chem/fingerprintStrategy/AtomPairsFngpr.o: chem/fingerprintStrategy/AtomPairsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/AtomPairsFngpr.o chem/fingerprintStrategy/AtomPairsFngpr.cpp

${OBJECTDIR}/chem/fingerprintStrategy/FingerprintStrategy.o: chem/fingerprintStrategy/FingerprintStrategy.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/FingerprintStrategy.o chem/fingerprintStrategy/FingerprintStrategy.cpp

${OBJECTDIR}/chem/fingerprintStrategy/MorganFngpr.o: chem/fingerprintStrategy/MorganFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/MorganFngpr.o chem/fingerprintStrategy/MorganFngpr.cpp

${OBJECTDIR}/chem/fingerprintStrategy/TopolLayeredFngpr1.o: chem/fingerprintStrategy/TopolLayeredFngpr1.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/TopolLayeredFngpr1.o chem/fingerprintStrategy/TopolLayeredFngpr1.cpp

${OBJECTDIR}/chem/fingerprintStrategy/TopolLayeredFngpr2.o: chem/fingerprintStrategy/TopolLayeredFngpr2.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/TopolLayeredFngpr2.o chem/fingerprintStrategy/TopolLayeredFngpr2.cpp

${OBJECTDIR}/chem/fingerprintStrategy/TopolSingleFngpr.o: chem/fingerprintStrategy/TopolSingleFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/TopolSingleFngpr.o chem/fingerprintStrategy/TopolSingleFngpr.cpp

${OBJECTDIR}/chem/fingerprintStrategy/TopolTorsFngpr.o: chem/fingerprintStrategy/TopolTorsFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/TopolTorsFngpr.o chem/fingerprintStrategy/TopolTorsFngpr.cpp

${OBJECTDIR}/chem/fingerprintStrategy/VectorFpFngpr.o: chem/fingerprintStrategy/VectorFpFngpr.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/fingerprintStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/fingerprintStrategy/VectorFpFngpr.o chem/fingerprintStrategy/VectorFpFngpr.cpp

${OBJECTDIR}/chem/morphing/Morphing.o: chem/morphing/Morphing.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphing/Morphing.o chem/morphing/Morphing.cpp

${OBJECTDIR}/chem/morphing/MorphingData.o: chem/morphing/MorphingData.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphing/MorphingData.o chem/morphing/MorphingData.cpp

${OBJECTDIR}/chem/morphing/MorphingFtors.o: chem/morphing/MorphingFtors.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphing
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphing/MorphingFtors.o chem/morphing/MorphingFtors.cpp

${OBJECTDIR}/chem/morphingStrategy/OpAddAtom.o: chem/morphingStrategy/OpAddAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpAddAtom.o chem/morphingStrategy/OpAddAtom.cpp

${OBJECTDIR}/chem/morphingStrategy/OpAddBond.o: chem/morphingStrategy/OpAddBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpAddBond.o chem/morphingStrategy/OpAddBond.cpp

${OBJECTDIR}/chem/morphingStrategy/OpBondContraction.o: chem/morphingStrategy/OpBondContraction.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpBondContraction.o chem/morphingStrategy/OpBondContraction.cpp

${OBJECTDIR}/chem/morphingStrategy/OpBondReroute.o: chem/morphingStrategy/OpBondReroute.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpBondReroute.o chem/morphingStrategy/OpBondReroute.cpp

${OBJECTDIR}/chem/morphingStrategy/OpInterlayAtom.o: chem/morphingStrategy/OpInterlayAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpInterlayAtom.o chem/morphingStrategy/OpInterlayAtom.cpp

${OBJECTDIR}/chem/morphingStrategy/OpMutateAtom.o: chem/morphingStrategy/OpMutateAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpMutateAtom.o chem/morphingStrategy/OpMutateAtom.cpp

${OBJECTDIR}/chem/morphingStrategy/OpRemoveAtom.o: chem/morphingStrategy/OpRemoveAtom.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpRemoveAtom.o chem/morphingStrategy/OpRemoveAtom.cpp

${OBJECTDIR}/chem/morphingStrategy/OpRemoveBond.o: chem/morphingStrategy/OpRemoveBond.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/morphingStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/morphingStrategy/OpRemoveBond.o chem/morphingStrategy/OpRemoveBond.cpp

${OBJECTDIR}/chem/simCoefStrategy/AllBitSimCoef.o: chem/simCoefStrategy/AllBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/AllBitSimCoef.o chem/simCoefStrategy/AllBitSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/AsymmetricSimCoef.o: chem/simCoefStrategy/AsymmetricSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/AsymmetricSimCoef.o chem/simCoefStrategy/AsymmetricSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/BraunBlanquetSimCoef.o: chem/simCoefStrategy/BraunBlanquetSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/BraunBlanquetSimCoef.o chem/simCoefStrategy/BraunBlanquetSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/CosineSimCoef.o: chem/simCoefStrategy/CosineSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/CosineSimCoef.o chem/simCoefStrategy/CosineSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/DiceSimCoef.o: chem/simCoefStrategy/DiceSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/DiceSimCoef.o chem/simCoefStrategy/DiceSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/KulczynskiSimCoef.o: chem/simCoefStrategy/KulczynskiSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/KulczynskiSimCoef.o chem/simCoefStrategy/KulczynskiSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/McConnaugheySimCoef.o: chem/simCoefStrategy/McConnaugheySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/McConnaugheySimCoef.o chem/simCoefStrategy/McConnaugheySimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/OnBitSimCoef.o: chem/simCoefStrategy/OnBitSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/OnBitSimCoef.o chem/simCoefStrategy/OnBitSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/RusselSimCoef.o: chem/simCoefStrategy/RusselSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/RusselSimCoef.o chem/simCoefStrategy/RusselSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/SokalSimCoef.o: chem/simCoefStrategy/SokalSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/SokalSimCoef.o chem/simCoefStrategy/SokalSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/TanimotoSimCoef.o: chem/simCoefStrategy/TanimotoSimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/TanimotoSimCoef.o chem/simCoefStrategy/TanimotoSimCoef.cpp

${OBJECTDIR}/chem/simCoefStrategy/TverskySimCoef.o: chem/simCoefStrategy/TverskySimCoef.cpp 
	${MKDIR} -p ${OBJECTDIR}/chem/simCoefStrategy
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/chem/simCoefStrategy/TverskySimCoef.o chem/simCoefStrategy/TverskySimCoef.cpp

${OBJECTDIR}/coord/KamadaKawaiReducer.o: coord/KamadaKawaiReducer.cpp 
	${MKDIR} -p ${OBJECTDIR}/coord
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/coord/KamadaKawaiReducer.o coord/KamadaKawaiReducer.cpp

${OBJECTDIR}/coord/PcaReducer.o: coord/PcaReducer.cpp 
	${MKDIR} -p ${OBJECTDIR}/coord
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/coord/PcaReducer.o coord/PcaReducer.cpp

${OBJECTDIR}/coord/ReducerFactory.o: coord/ReducerFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/coord
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/coord/ReducerFactory.o coord/ReducerFactory.cpp

${OBJECTDIR}/core/JobManager.o: core/JobManager.cpp 
	${MKDIR} -p ${OBJECTDIR}/core
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/core/JobManager.o core/JobManager.cpp

${OBJECTDIR}/core/NeighborhoodGenerator.o: core/NeighborhoodGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/core
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/core/NeighborhoodGenerator.o core/NeighborhoodGenerator.cpp

${OBJECTDIR}/core/NeighborhoodTaskQueue.o: core/NeighborhoodTaskQueue.cpp 
	${MKDIR} -p ${OBJECTDIR}/core
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/core/NeighborhoodTaskQueue.o core/NeighborhoodTaskQueue.cpp

${OBJECTDIR}/core/PathFinder.o: core/PathFinder.cpp 
	${MKDIR} -p ${OBJECTDIR}/core
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/core/PathFinder.o core/PathFinder.cpp

${OBJECTDIR}/core/PathFinderContext.o: core/PathFinderContext.cpp 
	${MKDIR} -p ${OBJECTDIR}/core
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/core/PathFinderContext.o core/PathFinderContext.cpp

${OBJECTDIR}/extensions/SAScore.o: extensions/SAScore.cpp 
	${MKDIR} -p ${OBJECTDIR}/extensions
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/extensions/SAScore.o extensions/SAScore.cpp

${OBJECTDIR}/tests/MorphingTest.o: tests/MorphingTest.cpp 
	${MKDIR} -p ${OBJECTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DNETBEANS_HACK -DRCF_MULTI_THREADED -DRCF_NO_AUTO_INIT_DEINIT -DRCF_USE_BOOST_ASIO -DRCF_USE_BOOST_SERIALIZATION -DRCF_USE_BOOST_THREADS -DRCF_USE_ZLIB -DTBB_USE_DEBUG=0 -I./. -I../common -I../dependencies/boost -I../dependencies/rcf/include -I../dependencies/zlib -I../dependencies/rdkit/Code -I../dependencies/tbb/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tests/MorphingTest.o tests/MorphingTest.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libbackend.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
