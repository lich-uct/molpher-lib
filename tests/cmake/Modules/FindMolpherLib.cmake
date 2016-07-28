FIND_PATH(MOLPHER_INCLUDE_DIR data_structs/ExplorationTree.hpp
        ${MOLPHER_INSTALL_DIR}/include/molpher
        )
FIND_LIBRARY(MOLPHER_LIBRARY molpher
        ${MOLPHER_INSTALL_DIR}/lib
        )