find_package(Boost REQUIRED system serialization)

set(GRAPH_SOURCE_FILES
        ../Core/Graph.hpp
        ../AOG_Base/AOG_Edge.h
        ../AOG_Base/AOG_Vertex.h
        ../AOG_Base/Symbolic_Rule.h
        ../AOG_Base/AOG.h
        ../AOG_Base/Symbolic_State.h)

set(GRAPH_TEST_FILE_1
        grammar_example.cpp)

set(GRAPH_TEST_FILE_2
        1d_sentence_exp.cpp)

set(GRAPH_TEST_FILE_3
        KL_distance.cpp)

set(
        TEST_AOG_MERGE_SRCS
        test_aog_merge.cpp
)

set(
        TEST_AOG_SIMPLIFY_SRCS
        test_aog_simplify.cpp
)

include_directories(include ../AOG_Base/)

add_executable(grammar_exp ${GRAPH_SOURCE_FILES} ${GRAPH_TEST_FILE_1})
add_executable(kl_distance ${GRAPH_SOURCE_FILES} ${GRAPH_TEST_FILE_3})

target_link_libraries(grammar_exp ${Boost_LIBRARIES})
target_link_libraries(kl_distance ${Boost_LIBRARIES})