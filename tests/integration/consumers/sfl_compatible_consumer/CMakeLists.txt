cmake_minimum_required(VERSION 3.20)

project(MatCalSflCompatibleConsumer LANGUAGES CXX)

find_package(MatCal CONFIG REQUIRED)

add_executable(matcal_sfl_compatible_consumer main.cpp)
target_compile_features(matcal_sfl_compatible_consumer PRIVATE cxx_std_20)
target_link_libraries(matcal_sfl_compatible_consumer PRIVATE MatCal::Linalg)
