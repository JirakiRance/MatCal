cmake_minimum_required(VERSION 3.20)

project(MatCalAddSubdirectoryConsumer LANGUAGES CXX)

set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
add_subdirectory(../.. matcal-build)

add_executable(matcal_add_subdirectory_consumer main.cpp)
target_link_libraries(matcal_add_subdirectory_consumer
    PRIVATE
        MatCal::Roots
        MatCal::Interpolation
        MatCal::Polynomial
        MatCal::Linalg
        MatCal::Legacy
)
