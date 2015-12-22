PREFIX=""

# prepare directories
if [[ -z "$VAR" ]]; then
    PREFIX=`pwd`
else
    mkdir -p $PREFIX
    PREFIX=$(realpath $PREFIX)
    cp -r . $PREFIX
fi

function build_deps()
{
    echo "Molpher Install: Building boost..."
    cd $PREFIX/dependencies/
    ./build-deps-linux-backend-ia64.sh --boost
    echo "Molpher Install: Building rdkit..."
    ./build-deps-linux-backend-ia64.sh --rdkit
    echo "Molpher Install: Building zlib..."
    ./build-deps-linux-backend-ia64.sh --zlib
    echo "Molpher Install: Building rcf..."
    ./build-deps-linux-backend-ia64.sh --rcf
    echo "Molpher Install: Building tbb..."
    ./build-deps-linux-backend-ia64.sh --tbb
}

function build_molpher()
{
    echo "Molpher Install: Building Molpher library file..."
    cd $PREFIX/backend/
    make CONF=Linux64_Debug_library
}

#function build_molpher_lib()
#{
#    echo "Molpher Install: Building Molpher-lib..."
#    cd $PREFIX/molpher-lib/
#    make CONF=Linux64_Debug
#}

build_deps && build_molpher && echo "Molpher was installed successfully in $PREFIX!" && exit
echo "Molpher installation failed!"
