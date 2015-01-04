#!/bin/sh

# Prepare directories
mkdir ./molph-redist-win32
mkdir ./molph-redist-win32/molpher

# Prepare binaries
cp ./backend/dist/Release/MinGW-Windows/backend.exe ./molph-redist-win32/molpher/molpher-srv.exe
cp ./frontend/dist/Release/MinGW-Windows/frontend.exe ./molph-redist-win32/molpher/molpher-gui.exe
cp ./dependencies/tbb/build/windows_ia32_gcc_mingw_release/tbb.dll ./molph-redist-win32/molpher
cp ./dependencies/mingw/libgcc_s_dw2-1.dll ./molph-redist-win32/molpher
cp ./dependencies/mingw/libstdc++-6.dll ./molph-redist-win32/molpher
cp ./dependencies/indigo/indigo-depict.exe ./molph-redist-win32/molpher
upx --best ./molph-redist-win32/molpher/molpher-srv.exe
upx --best ./molph-redist-win32/molpher/molpher-gui.exe

# Prepare testing files (TODO this should be probably removed in future)
mkdir ./molph-redist-win32/molpher/Inputs
cp ./backend/TestFiles/CID_10635-CID_15951529.sdf ./molph-redist-win32/molpher/Inputs

# Prepare auxiliary documents
cp ./license.txt ./molph-redist-win32/molpher

# Prepare redistributable package
cd ./molph-redist-win32
zip -9 -r molph-redist-win32 ./*
cd ..
mv ./molph-redist-win32/molph-redist-win32.zip ./

# Cleanup
rm -rf ./molph-redist-win32
