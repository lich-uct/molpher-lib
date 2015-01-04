Linux intel64 backend

1) download dependencies from https://www.assembla.com/spaces/molpher/documents/tag/dev file 'molph-dev-env-linux-intel64.tar'

2) unpack the downloaded file into ./dependencies

3) use ./depedencies/build-deps-win-ia64.sh to build
	-boost
	-zlib
	-tbb
	-rdkit
	-rfc

4) build the backend by executing 'make CONF=Linux64_Release' in backend directory
	
5) run './build-redist-linux-ia64.sh' to create dist directory with Molpher backend files

6) use './dist/run.sh' to run Molpher backend

