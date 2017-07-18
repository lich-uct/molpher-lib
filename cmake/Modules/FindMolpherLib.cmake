FIND_PATH(MOLPHER_INCLUDE_DIR data_structs/ExplorationTree.hpp
        ${CMAKE_INCLUDE_PATH}/molpher
        )
FIND_LIBRARY(MOLPHER_LIBRARY molpher
        ${CMAKE_LIBRARY_PATH}
        )