#!/bin/sh

# Copyright (c) 2013 Petr Škoda
# Copyright (c) 2016 Martin Šícho
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

set -e

# -------------------------------------- COMMENT -------------------------------------- #

# Dependency links
TBB_LINK=http://www.threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb42_20130725oss_lin.tgz
BOOST_LINK=https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz/download
RDKIT_LINK=https://github.com/rdkit/rdkit/archive/Release_2017_09_1.tar.gz

# Dependency compilation flags
want_help=0
TBB=
BOOST=
RDKIT=

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
            RDKIT=$1
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
                boost)    BOOST=$2;;
                rdkit)    RDKIT=$2;;
            esac
        fi
    fi
    return 1
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

print_all()
{
    echo "The following dependencies will be built:"
    if [ $TBB -eq 1 ]
    then
        echo "TBB"
    fi
    if [ $BOOST -eq 1 ]
    then
        echo "Boost"
    fi
    if [ $RDKIT -eq 1 ]
    then
        echo "RDKit"
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

    # cd tbb
    # just use make - it will take care of all
    # make compiler=gcc arc=ia64 runtime=gcc

    # rm -r -f ./build/windows_ia32_gcc_mingw_debug
    # rm -f ./build/windows_ia32_gcc_mingw_release/*.o
    # rm -f ./build/windows_ia32_gcc_mingw_release/*.d
    # cd ..
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
    count=`ls -1 boost_1_*.tar.gz 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv boost_1_*.tar.gz boost.tar.gz
    fi

    # download if doesn't exist
    if [ ! -f ./boost.tar.gz ]
    then
        wget --output-document=boost.tar.gz $BOOST_LINK
    fi

    echo "Unpacking..."
    tar --extract --gzip --file=boost.tar.gz
    mv boost_1_* boost
    cd boost

    chmod +x bootstrap.sh

    # build boost
    ./bootstrap.sh gcc --without-icu
    ./b2 variant=release link=static,shared runtime-link=shared threading=multi toolset=gcc cxxflags=\"-fPIC\" --without-chrono --without-exception --without-graph --without-graph_parallel --without-iostreams --without-locale --without-math --without-mpi --without-python --without-random --without-test --without-timer --without-wave

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
    count=`ls -1 RDKit_20*.tar.gz 2>/dev/null | wc -l`
    if [ $count -ge 1 ]
    then
        mv RDKit_20*.tar.gz rdkit.tar.gz
    fi

    # download if doesn't exist
    if [ ! -f ./rdkit.tar.gz ]
    then
        wget --output-document=rdkit.tar.gz $RDKIT_LINK
    fi

    echo "Unpacking..."
    tar --extract --gzip --file=rdkit.tar.gz
    mv rdkit-Release_20* rdkit
    cd rdkit
    cmd=cmd
    # prepare script to build RDKit
    echo "RDBASE=`pwd`" > $cmd
    echo "LD_LIBRARY_PATH=`pwd`/lib:`pwd`/../boost/stage/lib" >> $cmd
    echo "mkdir build" >> $cmd
    echo "cd build" >> $cmd
    echo "cmake -G \"Unix Makefiles\" -D RDK_BUILD_PYTHON_WRAPPERS= -D BOOST_ROOT=../boost -D RDK_INSTALL_STATIC_LIBS=ON -D Boost_USE_STATIC_LIBS=ON .." >> $cmd
    echo "make -j 4" >> $cmd
    echo "make install" >> $cmd
    echo "cd .." >> $cmd
    echo "rm -r -f ./build" >> $cmd

    # executer script and remove file with script
    chmod +x $cmd
    ./$cmd
    rm $cmd

    cd ..
}

setval_all 0

if [ $# -eq 0 ]
then
    echo "To show help use './build_deps.sh -h'"
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
        -tbb | --tbb)
            TBB=1
        ;;
        -boost | --boost)
            BOOST=1
        ;;
        -rdkit | --rdkit)
            RDKIT=1
        ;;
    esac
done

if [ $want_help -eq 1 ]
then
    echo "Options:"
    echo "-help | --help | -h                  this help"
    echo "-all | --all | -a                    build all dependencies"
    echo "-tbb | --tbb                         build TBB"
    echo "-boost | --boost                     build Boost"
    echo "-rdkit | --rdkit                     build RDKit"
fi

# Call build functions

rdkit_dep

print_all

build_tbb
build_boost
build_rdkit

exit
