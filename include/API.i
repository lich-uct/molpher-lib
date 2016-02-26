%include <std_shared_ptr.i>;
%include "std_unique_ptr.i"

%shared_ptr(ExplorationTree);
%shared_ptr(MolpherMol);

// tree operations
%include "operations/operations.i"

// callbacks
//%include "operations/callbacks/callbacks.i"

// data structs
%include "data_structs/data_structs.i"