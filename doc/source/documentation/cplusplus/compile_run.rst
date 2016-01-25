Compiling and running your programs
===================================

Unfortunately, because the current API implementation relies heavily on the original Molpher code,
additional header files from the `Molpher-backend project`, its dependencies and
some `common headers` need to be included during compilation.
Therefore, if you are using the headers located in :file:`molpher-lib/include/molpher_API`,
make sure to tell the compiler to include the following directories (all paths are relative to the repository root):

   - :file:`molpher-lib/include/`
   - :file:`backend/`
   - :file:`common/`
   - :file:`dependencies/rdkit/Code`
   - :file:`dependencies/boost`
   - :file:`dependencies/tbb/include`

Also note that the library is dynamically linked against your program during runtime.
Therefore, it is necessary to include the :file:`molpher-lib/lib/` directory
with the compiled library files. You can use the :command:`LD_LIBRARY_PATH`
environment variable for that.