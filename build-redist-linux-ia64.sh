#!/bin/sh

# delete directory
rm -f -r dist

#create directory
mkdir dist

# copy molpher-srv and SASrore.dat
cp backend/dist/*/*/molpher-srv dist/molpher-srv
cp backend/SAScore.dat dist/SAScore.dat

# copy libraries
mkdir dist/libs/
cp dependencies/tbb/lib/intel64/gcc4.4/libtbb.so.2 dist/libs/libtbb.so.2
cp dependencies/tbb/lib/intel64/gcc4.4/libtbbmalloc.so.2 dist/libs/libtbbmalloc.so.2
cp dependencies/tbb/lib/intel64/gcc4.4/libtbbmalloc_proxy.so.2 dist/libs/libtbbmalloc_proxy.so.2

# create run file
echo '#!/bin/sh' > dist/run.sh
echo 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(pwd)/libs/' >> dist/run.sh
echo './molpher-srv' >> dist/run.sh

chmod +x dist/run.sh

