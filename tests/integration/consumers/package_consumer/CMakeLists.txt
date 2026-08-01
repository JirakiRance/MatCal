cmake_minimum_required(VERSION 3.20)

project(MatCalPackageConsumer LANGUAGES CXX)

find_package(MatCal 0.3 CONFIG REQUIRED)

foreach(matcal_target IN ITEMS
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
    if(NOT TARGET ${matcal_target})
        message(FATAL_ERROR "Missing exported MatCal target: ${matcal_target}")
    endif()
endforeach()

get_target_property(matcal_legacy_links MatCal::Legacy INTERFACE_LINK_LIBRARIES)
foreach(expected_dependency IN ITEMS
    MatCal::Polynomial
    MatCal::Roots
    MatCal::Interpolation
    MatCal::Calculus
    MatCal::ODE
    MatCal::Nonlinear
    MatCal::LeastSquares
)
    if(NOT matcal_legacy_links MATCHES "${expected_dependency}")
        message(FATAL_ERROR "MatCal::Legacy does not export dependency ${expected_dependency}")
    endif()
endforeach()

get_target_property(matcal_linalg_links MatCal::Linalg INTERFACE_LINK_LIBRARIES)
if(matcal_linalg_links MATCHES "Legacy")
    message(FATAL_ERROR "MatCal::Linalg must not depend on MatCal::Legacy")
endif()

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
