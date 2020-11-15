find_package(OpenCV REQUIRED )
find_package(Boost REQUIRED system serialization)

set(
    SRCS_Sample
    src/components.cpp
    src/sample.cpp
    src/utils.cpp
    src/face_params.cpp
)

set(
    SRCS_Learn
    src/components.cpp
    src/learn.cpp
    src/utils.cpp
    src/face_params_learn.cpp
)


include_directories(
    ${OpenCV_INCLUDE_DIRS}
    include
    ../AOG_Base/
    ../Learner/Parameter_Learner/
)

MESSAGE(STATUS "Include dirs ${OpenCV_INCLUDE_DIRS}")
MESSAGE(STATUS "LINK LIBRARIES ${OpenCV_LIBS}")

add_executable(cartoon_face ${SRCS_Sample})
add_executable(parameter_learning ${SRCS_Learn})

target_link_libraries(cartoon_face ${Boost_LIBRARIES} ${OpenCV_LIBS})
target_link_libraries(parameter_learning ${Boost_LIBRARIES} ${OpenCV_LIBS})
