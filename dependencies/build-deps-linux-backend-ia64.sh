#!/bin/sh

# -------------------------------------- COMMENT -------------------------------------- #

# Dependency links
TBB_LINK=http://www.threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb42_20130725oss_lin.tgz
BOOST_LINK=http://sourceforge.net/projects/boost/files/boost/1.49.0/boost_1_49_0.tar.gz/download
RDKIT_LINK=https://github.com/rdkit/rdkit/archive/Release_2014_03_1.tar.gz

# version 1.2.7 is not aviable anymore
ZLIB_LINK=http://zlib.net/zlib-1.2.7.tar.gz
RCF_LINK=http://www.deltavsoft.com/downloads/RCF-1.3.1.tar.gz

# Dependency compilation flags
want_help=0
TBB=
BOOST=
BOOSTMPI=0
RDKIT=
MPICH2=0
ZLIB=
RCF=

# Toggle compilation of all dependencies.
# $1   compilation flag (1 - build all, 0 - don't build anything)
setval_all()
{
    if [ $# -eq 1 ]
    then
        if [ $1 -eq 0 -o $1 -eq 1 ]
        then
            TBB=$1
            BOOST=$1
            #BOOSTMPI=$1
            RDKIT=$1
            #MPICH2=$1
            ZLIB=$1
            RCF=$1
            return 0
        fi
    fi
    return 1
}

# Toggle compilation of a particular dependency.
# $1   dependency name
# $2   compilation flag (1 - build, 0 - don't build)
setval()
{
    if [ $# -eq 2 ]
    then
        if [ $2 -eq 0 -o $2 -eq 1 ]
        then
            case $1 in
                tbb)      TBB=$2;;
                qt)       QT=$2;;
                boost)    BOOST=$2;;
                boostmpi)
                          BOOST=$2
                          BOOSTMPI=$2
                          MPICH2=$2
                ;;
                rdkit)    RDKIT=$2;;
                mpich2)   MPICH2=$2;;
                zlib)     ZLIB=$2;;
                rcf)      RCF=$2;;
                mingw)    MINGW=$2;;
                indigo)   INDIGO=$2;;
            esac
        fi
    fi
    return 1
}

# Boost with MPI support is dependent on MPICH2.
boostmpi_dep()
{
    if [ $BOOSTMPI -eq 1 ]
    then
        if [ ! -d ./mpich2/lib -a $MPICH2 -eq 0 ]
        then
            echo "Boost MPI is dependent on MPICH2 libraries, which were not found."
            while [ 1 -eq 1 ]
            do
                echo "Do you want to build MPICH2? (y/n)"
                read resp
                resp=${resp^^}
                case $resp in
                    Y)
                        MPICH2=1
                        break
                    ;;
                    N)
                        echo "Boost MPI will not be built."
                        BOOSTMPI=0
                        break
                    ;;
                    *)
                        echo "Bad response"
                    ;;
                esac
            done
        fi
    fi
}

# MPICH2 must be already installed on the system.
mpich2_dep()
{
    if [ $MPICH2 -eq 1 ]
    then
        if [ -z "$MPICH2_PATH" -o ! -d "$MPICH2_PATH/lib" ]
        then
            echo "MPICH2 installation not found."
            MPICH2=0
            if [ $BOOSTMPI -eq 1 ]
            then
                echo "Boost MPI will not be built."
                BOOSTMPI=0
            fi
        fi
    fi
}

# RDKit is dependent on Boost.
rdkit_dep()
{
    if [ $RDKIT -eq 1 ]
    then
        if [ ! -d ./boost/stage/lib -a $BOOST -eq 0 ] 
        then
            echo "RDKit is dependent on Boost libraries, which were not found."
            while [ 1 -eq 1 ]
            do
                echo "Do you want to build Boost? (y/n)"
                read resp
                resp=${resp^^}
                case $resp in
                    Y)
                        BOOST=1
                        break
                    ;;
                    N)
                        echo "RDKit will not be built."
                        RDKIT=0
                        break
                    ;;
                    *)
                        echo "Bad response"
                    ;;
                esac
            done
        fi
    fi
}

# RCF is dependent on Boost and zlib.
rcf_dep()
{
    if [ $RCF -eq 1 ]
    then
        if [ ! -d ./boost/stage/lib -a $BOOST -eq 0 ]
        then
            echo "RCF is dependent on Boost libraries, which were not found."
            while [ 1 -eq 1 ]
            do
                echo "Do you want to build Boost? (y/n)"
                read resp
                resp=${resp^^}
                case $resp in
                    Y)
                        BOOST=1
                        break
                    ;;
                    N)
                        echo "RCF will not be built."
                        RCF=0
                        break
                    ;;
                    *)
                        echo "Bad response"
                    ;;
                esac
            done
        fi
        if [ ! -e ./zlib/libz.a -a $ZLIB -eq 0 ]
        then
            echo "RCF is dependent on zlib library, which was not found."
            while [ 1 -eq 1 ]
            do
                echo "Do you want to build zlib? (y/n)"
                read resp
                resp=${resp^^}
                case $resp in
                    Y)
                        ZLIB=1
                        break
                    ;;
                    N)
                        echo "RCF will not be built"
                        RCF=0
                        break
                    ;;
                    *)
                        echo "Bad response"
                    ;;
                esac
            done
        fi
    fi
}

print_all()
{
    echo "The following dependencies will be built:"
    if [ $TBB -eq 1 ] 
    then
        echo "TBB"
    fi
    if [ $BOOST -eq 1 ]
    then
        if [ $BOOSTMPI -eq 1 ]
        then
            echo "Boost with MPI"
        else
            echo "Boost"
        fi
    fi
    if [ $RDKIT -eq 1 ]
    then
        echo "RDKit"
    fi
    if [ $MPICH2 -eq 1 ]
    then
        echo "MPICH2"
    fi
    if [ $ZLIB -eq 1 ]
    then
        echo "zlib"
    fi
    if [ $RCF -eq 1 ]
    then
        echo "RCF" 
    fi
}

build_tbb()
{
    if [ $TBB -eq 0 ]
    then
        return
    fi
    
    echo "Building TBB"
    
    if [ -d ./tbb ]
    then 
        echo "Cleaning..."
        rm -r -f ./tbb
    fi
    
    # changing the name of the package
    count=`ls -1 tbb42_*.tgz 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv tbb42_*.tgz tbb.tar.gz
    fi
    
    # download if doesn't exist
    if [ ! -f ./tbb.tar.gz ]
    then
        wget --output-document=tbb.tar.gz $TBB_LINK
    fi
    
    echo "Unpacking..."
    tar --extract --gzip --file=tbb.tar.gz
    mv tbb42_*oss tbb

    #cd tbb
    # just use make - it will take care of all
    #make compiler=gcc arc=ia64 runtime=gcc
    
    ##rm -r -f ./build/windows_ia32_gcc_mingw_debug
    ##rm -f ./build/windows_ia32_gcc_mingw_release/*.o
    ##rm -f ./build/windows_ia32_gcc_mingw_release/*.d
    cd ..
}

build_boost()
{
    if [ $BOOST -eq 0 ]
    then
        return
    fi
    
    echo "Building Boost"
    #return
    
    if [ -d ./boost ] 
    then 
        echo "Cleaning..."
        rm -r -f ./boost
    fi
    
    # changing the name of the package
    count=`ls -1 boost_1_4*.tar.gz 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv boost_1_4*.tar.gz boost.tar.gz
    fi
    
    # download if doesn't exist
    if [ ! -f ./boost.tar.gz ]
    then
        wget --output-document=boost.tar.gz $BOOST_LINK
    fi
    
    echo "Unpacking..."
    tar --extract --gzip --file=boost.tar.gz
    mv boost_1_4* boost
    cd boost
    
    chmod +x bootstrap.sh
    
    if [ $BOOSTMPI -eq 1 ]
    then
        echo "using mpi : : <include>../mpich2/include <library-path>../mpich2/lib <find-static-library>cxx <find-static-library>mpi ;" >> ./tools/build/v2/user-config.jam
    fi

    # build boost 
    if [ $BOOSTMPI -eq 1 ]
    then
        ./bootstrap.sh gcc --without-icu
        ./b2 variant=release link=static runtime-link=shared threading=multi toolset=gcc cxxflags=\"-fpermissive -Wno-unused-value -fPIC\" --without-chrono --without-exception --without-graph --without-graph_parallel --without-iostreams --without-locale --without-math --without-python --without-random --without-test --without-timer --without-wave
    else
        ./bootstrap.sh gcc --without-icu
        ./b2 variant=release link=static runtime-link=shared threading=multi toolset=gcc cxxflags=\"-fPIC\" --without-chrono --without-exception --without-graph --without-graph_parallel --without-iostreams --without-locale --without-math --without-mpi --without-python --without-random --without-test --without-timer --without-wave
    fi
    
    rm -r -f ./bin.v2
    cd ..
}

build_rdkit()
{
    if [ $RDKIT -eq 0 ]
    then
        return
    fi
    
    echo "Building RDKit"
    
    if [ -d ./rdkit ] 
    then 
        echo "Cleaning..."
        rm -r -f ./rdkit
    fi
    
    # changing the name of the package
    count=`ls -1 rdkit-*.* 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv rdkit-*.* rdkit.tar.gz
    fi
    
    # download if doesn't exist
    if [ ! -f ./rdkit.tar.gz ]
    then
        wget --output-document=rdkit.tar.gz $RDKIT_LINK
    fi
    
    echo "Unpacking..."
    tar --extract --gzip --file=rdkit.tar.gz
    mv rdkit-* rdkit
    cd rdkit
    cmd=cmd
    # prepare script to build RDKit
    echo "RDBASE=`pwd`" > $cmd
    echo "LD_LIBRARY_PATH=`pwd`/lib:`pwd`/../boost/stage/lib" >> $cmd
    echo "mkdir build" >> $cmd
    echo "cd build" >> $cmd
    echo "cmake -G \"Unix Makefiles\" -D RDK_BUILD_PYTHON_WRAPPERS= -D BOOST_ROOT=../boost -D RDK_INSTALL_STATIC_LIBS=ON -D Boost_USE_STATIC_LIBS=ON .." >> $cmd
    echo "make" >> $cmd
    echo "make install" >> $cmd
    echo "cd .." >> $cmd
    echo "rm -r -f ./build" >> $cmd

    # executer script and remove file with script
    chmod +x $cmd
    ./$cmd
    rm $cmd
    
    cd ..
}

build_mpich2()
{
    if [ $MPICH2 -eq 0 ]
    then
        return
    fi
    
    echo "Building MPICH2"
	echo "Require MPICH installation."
    
    if [ -d ./mpich2 ] 
    then 
        echo "Cleaning..."
        rm -r -f ./mpich2
    fi
    
    mkdir mpich2
    cp --recursive "$MPICH2_PATH/include" ./mpich2/include
    cp --recursive "$MPICH2_PATH/lib" ./mpich2/lib
}

build_zlib()
{
    if [ $ZLIB -eq 0 ]
    then
        return
    fi
    
    echo "Building zlib"
    #return
    
    if [ -d ./zlib ] 
    then 
        echo "Cleaning..."
        rm -r -f ./zlib
    fi
    
    # changing the name of the package
    count=`ls -1 zlib-1.2*.tar.gz 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv zlib-1.2*.tar.gz zlib.tar.gz
    fi
    
    # download if doesn't exist
    if [ ! -f ./zlib.tar.gz ]
    then
        wget --output-document=zlib.tar.gz $ZLIB_LINK
    fi
    
    echo "Unpacking..."
    tar --extract --gzip --file=zlib.tar.gz
    mv zlib-1.2* zlib
    cd zlib
    
    # script to build zlib
    cmd=cmd 
	echo "CFLAGS='-mstackrealign -fPIC -O3' ./configure --prefix=. --shared" > $cmd
	echo "make test" >> $cmd
	echo "make install" >> $cmd
	
	chmod +x $cmd
    ./$cmd
    rm $cmd
    
    cd ..
}

build_rcf()
{
    if [ $RCF -eq 0 ]
    then
        return
    fi
    
    echo "Building RCF"
    
    if [ -d ./rcf ] 
    then 
        echo "Cleaning..."
        rm -r -f ./rcf
    fi
    
    # changing the name of the package
    count=`ls -1 RCF-1.3*.tar.gz 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv RCF-1.3*.tar.gz rcf.tar.gz
    fi
    
    # download if doesn't exist
    if [ ! -f ./rcf.tar.gz ]
    then
        wget --output-document=rcf.tar.gz $RCF_LINK
    fi
    
    echo "Unpacking..."
    tar --extract --gzip --file=rcf.tar.gz
    mv RCF-1.3* rcf
    cd rcf
    
    pwd
    
    # RCF Makefile
    echo 'SET(CMAKE_C_COMPILER gcc-4.6)' > CMakeLists.txt
    echo 'SET(CMAKE_CXX_COMPILER g++-4.6)' >> CMakeLists.txt
    echo 'CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )' >> CMakeLists.txt
    echo 'PROJECT( RCF )' >> CMakeLists.txt
    echo 'ADD_DEFINITIONS( -DBOOST_ALL_NO_LIB -DBOOST_THREAD_USE_LIB -DRCF_MULTI_THREADED -DRCF_USE_BOOST_THREADS -DRCF_USE_BOOST_ASIO -DRCF_USE_ZLIB -DRCF_USE_BOOST_SERIALIZATION -DRCF_NO_AUTO_INIT_DEINIT -Wno-deprecated -Wno-attributes -Wno-write-strings -O2 )' >> CMakeLists.txt
    echo 'INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/../boost ${CMAKE_SOURCE_DIR}/../zlib ${CMAKE_SOURCE_DIR}/include )' >> CMakeLists.txt
    echo 'SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin )' >> CMakeLists.txt
    echo 'SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin )' >> CMakeLists.txt
    echo 'SET(CMAKE_BUILD_TYPE Release)' >> CMakeLists.txt
    echo 'SET(Boost_USE_STATIC_LIBS ON)' >> CMakeLists.txt
    echo 'SET(CMAKE_VERBOSE_MAKEFILE ON)' >> CMakeLists.txt
    echo 'SET(CMAKE_CXX_FLAGS "-fPIC")' >> CMakeLists.txt
    echo 'LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/bin )' >> CMakeLists.txt
    echo 'ADD_LIBRARY( RcfLib STATIC ${CMAKE_SOURCE_DIR}/src/RCF/RCF.cpp )' >> CMakeLists.txt
    
    # script to build RCF with MinGW
	cmd=cmd
    echo "cmake -G \"Unix Makefiles\" ." > $cmd
    echo "make" >> $cmd
    echo "strip --strip-unneeded ./bin/libRcfLib.a" >> $cmd
    chmod +x $cmd 
	./$cmd
    rm $cmd
    
    rm -r -f ./CMakeFiles
    
    cd ..
}

setval_all 0

if [ $# -eq 0 ]
then
    echo "To show help use 'build-deps-win-ia32.sh -h'"
    exit
fi

for option
do
    option=`echo $option | tr '[A-Z]' '[a-z]'`
    case $option in
        -help | --help | -h)
            want_help=1
        ;;
        -all | --all | -a)
            setval_all 1
        ;;
        -with=* | --with=*)
            setval_all 0
            build_list=`expr "x$option" : "x-*with=\(.*\)"`
            old_IFS=$IFS
            IFS=,
            for build in $build_list
            do
                setval $build 1
            done
            IFS=$old_IFS
        ;;
        -without=* | --without=*)
            setval_all 1
            build_list=`expr "x$option" : "x-*without=\(.*\)"`
            old_IFS=$IFS
            IFS=,
            for build in $build_list
            do
                setval $build 0
            done
            IFS=$old_IFS
        ;;
        -tbb | --tbb)
            TBB=1
        ;;
        -boost | --boost)
            BOOST=1
        ;;
        -boostmpi | --boostmpi)
            BOOST=1
            BOOSTMPI=1
        ;;
        -rdkit | --rdkit)
            RDKIT=1
        ;;
        -mpich2 | --mpich2)
            MPICH2=1
        ;;
        -zlib | --zlib)
            ZLIB=1
        ;;
        -rcf | --rcf)
            RCF=1
        ;;
    esac
done

if [ $want_help -eq 1 ]
then
    echo "Options:"
    echo "-help | --help | -h                  this help"
    echo "-all | --all | -a                    build all dependencies"
    echo "                                     except MPICH2 and Boost MPI"
    echo "-with=list | --with=list             build all dependencies in the list"
    echo "                                     dependencies are separated with ','"
    echo "-without=list | --without=list       build all dependencies "
    echo "                                     except MPICH2 and Boost MPI"
    echo "                                     without the dependencies in the list"
    echo "                                     dependencies are separated with ','"
    echo "-tbb | --tbb                         build TBB"
    echo "-boost | --boost                     build Boost (without MPI support)"
    echo "-boostmpi | --boostmpi               build Boost (with MPI support)"
    echo "-rdkit | --rdkit                     build RDKit"
    echo "-mpich2 | --mpich2                   build MPICH2"
    echo "-zlib | --zlib                       build zlib"
    echo "-rcf | --rcf                         build RCF"
    echo
    echo "Examples:"
    echo "'build-deps-win-ia32.sh -with=boost,rdkit' will build Boost and RDKit"
    echo "'build-deps-win-ia32.sh -all -boostmpi' will build all dependencies"
    echo "                                        including Boost MPI and MPICH2"
    echo
fi

# Add cmake - used for metacentrum
module add cmake-2.8

# Call build functions

boostmpi_dep
mpich2_dep
rdkit_dep
rcf_dep

print_all

build_tbb
build_mpich2
build_boost
build_rdkit
build_zlib
build_rcf

exit
