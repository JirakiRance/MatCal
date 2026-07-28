cmake_minimum_required(VERSION 3.20)

project(MatCalPackageConsumer LANGUAGES CXX)

find_package(MatCal CONFIG REQUIRED)

add_executable(matcal_package_consumer main.cpp)
target_compile_features(matcal_package_consumer PRIVATE cxx_std_20)
target_link_libraries(matcal_package_consumer PRIVATE
    MatCal::Legacy
    MatCal::Linalg
    MatCal::Polynomial
    MatCal::Roots
    MatCal::Interpolation
    MatCal::Calculus
    MatCal::ODE
    MatCal::Nonlinear
    MatCal::LeastSquares
)
