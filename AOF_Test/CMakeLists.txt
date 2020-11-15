find_package(Boost REQUIRED)

set(GRAPH_SOURCE_FILES
        ../Core/Graph.hpp
        ../AOG_Base/AOG_Edge.h
        ../AOG_Base/AOG_Vertex.h
        ../AOG_Base/Symbolic_Rule.h
        ../AOG_Base/AOG.h
        ../AOG_Base/Symbolic_State.h)

set(STRUCTURE_LEARNER
        test_learner.cpp)

include_directories(include
                    ../AOG_Base/
                    ../Learner/Structure_Learner
                    ../Utils/)

add_executable(aof_learner ${GRAPH_SOURCE_FILES} ${STRUCTURE_LEARNER})
add_executable(grammar_exp1 ${GRAPH_SOURCE_FILES} grammar_exp1.cpp)

target_link_libraries(aof_learner ${Boost_LIBRARIES})