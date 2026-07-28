cmake_minimum_required(VERSION 3.20)

project(MatCal
    VERSION 0.1.0
    DESCRIPTION "Legacy C++ numerical calculation library"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(MATCAL_ENABLE_WARNINGS "Enable compiler warnings for MatCal test targets" ON)
option(MATCAL_ENABLE_SANITIZERS "Enable AddressSanitizer/UndefinedBehaviorSanitizer for MatCal test targets when supported" OFF)

add_library(matcal_polynomial INTERFACE)
add_library(MatCal::Polynomial ALIAS matcal_polynomial)

target_compile_features(matcal_polynomial INTERFACE cxx_std_20)
target_include_directories(matcal_polynomial
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_polynomial PROPERTIES EXPORT_NAME Polynomial)

add_library(matcal_roots INTERFACE)
add_library(MatCal::Roots ALIAS matcal_roots)

target_compile_features(matcal_roots INTERFACE cxx_std_20)
target_include_directories(matcal_roots
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_roots PROPERTIES EXPORT_NAME Roots)

add_library(matcal_interpolation INTERFACE)
add_library(MatCal::Interpolation ALIAS matcal_interpolation)

target_compile_features(matcal_interpolation INTERFACE cxx_std_20)
target_include_directories(matcal_interpolation
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(matcal_interpolation INTERFACE MatCal::Polynomial)
set_target_properties(matcal_interpolation PROPERTIES EXPORT_NAME Interpolation)

add_library(matcal_calculus INTERFACE)
add_library(MatCal::Calculus ALIAS matcal_calculus)

target_compile_features(matcal_calculus INTERFACE cxx_std_20)
target_include_directories(matcal_calculus
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_calculus PROPERTIES EXPORT_NAME Calculus)

add_library(matcal_ode INTERFACE)
add_library(MatCal::ODE ALIAS matcal_ode)

target_compile_features(matcal_ode INTERFACE cxx_std_20)
target_include_directories(matcal_ode
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_ode PROPERTIES EXPORT_NAME ODE)

add_library(matcal_nonlinear INTERFACE)
add_library(MatCal::Nonlinear ALIAS matcal_nonlinear)

target_compile_features(matcal_nonlinear INTERFACE cxx_std_20)
target_include_directories(matcal_nonlinear
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_nonlinear PROPERTIES EXPORT_NAME Nonlinear)

add_library(matcal_least_squares INTERFACE)
add_library(MatCal::LeastSquares ALIAS matcal_least_squares)

target_compile_features(matcal_least_squares INTERFACE cxx_std_20)
target_include_directories(matcal_least_squares
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_least_squares PROPERTIES EXPORT_NAME LeastSquares)

add_library(matcal_legacy INTERFACE)
add_library(MatCal::Legacy ALIAS matcal_legacy)

target_compile_features(matcal_legacy INTERFACE cxx_std_20)
target_include_directories(matcal_legacy
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(matcal_legacy INTERFACE
    MatCal::Polynomial
    MatCal::Roots
    MatCal::Interpolation
    MatCal::Calculus
    MatCal::ODE
    MatCal::Nonlinear
    MatCal::LeastSquares
)
set_target_properties(matcal_legacy PROPERTIES EXPORT_NAME Legacy)

add_library(matcal_linalg
    src/linalg/DenseSolver.cpp
    src/linalg/SymmetricSkylineMatrix.cpp
    src/linalg/SkylineLdlt.cpp
)
add_library(MatCal::Linalg ALIAS matcal_linalg)

target_compile_features(matcal_linalg PUBLIC cxx_std_20)
target_include_directories(matcal_linalg
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_linalg PROPERTIES EXPORT_NAME Linalg)
target_link_libraries(matcal_nonlinear INTERFACE MatCal::Linalg)
target_link_libraries(matcal_least_squares INTERFACE MatCal::Linalg MatCal::Polynomial)

function(matcal_configure_test_target target_name)
    if(MATCAL_ENABLE_WARNINGS)
        if(MSVC)
            target_compile_options(${target_name} PRIVATE /W4 /permissive-)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
        endif()
    endif()

    if(MATCAL_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target_name} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target_name} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()

function(matcal_configure_linalg_target target_name)
    if(MATCAL_ENABLE_WARNINGS)
        if(MSVC)
            target_compile_options(${target_name} PRIVATE /W4 /permissive- /WX)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic -Werror)
        endif()
    endif()

    if(MATCAL_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target_name} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target_name} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()

matcal_configure_linalg_target(matcal_linalg)

install(TARGETS matcal_polynomial matcal_roots matcal_interpolation matcal_calculus matcal_ode matcal_nonlinear matcal_least_squares matcal_legacy matcal_linalg
    EXPORT MatCalTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES
    src/Basics.hpp
    src/Insert.hpp
    src/Iteration.hpp
    src/Matrix.hpp
    src/QinJiuShao.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT MatCalTargets
    NAMESPACE MatCal::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MatCal
)

configure_package_config_file(
    cmake/MatCalConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MatCalConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MatCal
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MatCalConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MatCalConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MatCalConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MatCal
)

include(CTest)

if(BUILD_TESTING)
    function(matcal_add_legacy_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    matcal_add_legacy_test(matcal_characterization_tests
        tests/characterization_tests.cpp
    )

    matcal_add_legacy_test(matcal_legacy_triangular_regression_tests
        tests/legacy_triangular_regression_tests.cpp
    )

    matcal_add_legacy_test(matcal_legacy_zero_size_regression_tests
        tests/legacy_zero_size_regression_tests.cpp
    )

    matcal_add_legacy_test(matcal_legacy_numerics_regression_tests
        tests/legacy_numerics_regression_tests.cpp
    )

    matcal_add_legacy_test(matcal_legacy_sor_regression_tests
        tests/legacy_sor_regression_tests.cpp
    )

    matcal_add_legacy_test(matcal_legacy_matrix_linalg_adapter_tests
        tests/legacy_matrix_linalg_adapter_tests.cpp
    )

    function(matcal_add_linalg_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Linalg)
        matcal_configure_linalg_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_polynomial_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Polynomial MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_roots_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Roots MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_interpolation_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Interpolation MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_calculus_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Calculus MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_ode_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::ODE MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_nonlinear_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::Nonlinear MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    function(matcal_add_least_squares_test target_name)
        add_executable(${target_name} ${ARGN})
        target_link_libraries(${target_name} PRIVATE MatCal::LeastSquares MatCal::Legacy)
        matcal_configure_test_target(${target_name})
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    matcal_add_polynomial_test(matcal_polynomial_tests
        tests/polynomial_tests.cpp
    )

    matcal_add_polynomial_test(matcal_polynomial_public_header_tests
        tests/polynomial_public_headers_main.cpp
        tests/header_polynomial_self_contained.cpp
    )

    matcal_add_polynomial_test(matcal_polynomial_multi_tu_tests
        tests/polynomial_multi_tu_main.cpp
        tests/polynomial_multi_tu_a.cpp
        tests/polynomial_multi_tu_b.cpp
    )

    matcal_add_roots_test(matcal_roots_tests
        tests/roots_tests.cpp
    )

    matcal_add_interpolation_test(matcal_interpolation_tests
        tests/interpolation_tests.cpp
    )

    matcal_add_interpolation_test(matcal_interpolation_remaining_tests
        tests/interpolation_remaining_tests.cpp
    )

    add_executable(matcal_roots_interpolation_public_header_tests
        tests/roots_interpolation_public_headers_main.cpp
        tests/header_roots_self_contained.cpp
        tests/header_linear_interpolator_self_contained.cpp
        tests/header_cubic_spline_self_contained.cpp
        tests/header_polynomial_interpolation_self_contained.cpp
    )
    target_link_libraries(matcal_roots_interpolation_public_header_tests PRIVATE
        MatCal::Roots
        MatCal::Interpolation
    )
    matcal_configure_test_target(matcal_roots_interpolation_public_header_tests)
    add_test(NAME matcal_roots_interpolation_public_header_tests COMMAND matcal_roots_interpolation_public_header_tests)

    add_executable(matcal_roots_interpolation_multi_tu_tests
        tests/roots_interpolation_multi_tu_main.cpp
        tests/roots_multi_tu_a.cpp
        tests/interpolation_multi_tu_b.cpp
        tests/interpolation_remaining_multi_tu_c.cpp
    )
    target_link_libraries(matcal_roots_interpolation_multi_tu_tests PRIVATE
        MatCal::Roots
        MatCal::Interpolation
    )
    matcal_configure_test_target(matcal_roots_interpolation_multi_tu_tests)
    add_test(NAME matcal_roots_interpolation_multi_tu_tests COMMAND matcal_roots_interpolation_multi_tu_tests)

    add_executable(matcal_pt_style_consumer_tests
        tests/pt_style_consumer/main.cpp
    )
    target_link_libraries(matcal_pt_style_consumer_tests PRIVATE MatCal::Legacy)
    matcal_configure_test_target(matcal_pt_style_consumer_tests)
    add_test(NAME matcal_pt_style_consumer_tests COMMAND matcal_pt_style_consumer_tests)

    matcal_add_calculus_test(matcal_calculus_tests
        tests/calculus_tests.cpp
    )

    matcal_add_ode_test(matcal_ode_tests
        tests/ode_tests.cpp
    )

    add_executable(matcal_calculus_ode_public_header_tests
        tests/calculus_ode_public_headers_main.cpp
        tests/header_calculus_self_contained.cpp
        tests/header_ode_self_contained.cpp
    )
    target_link_libraries(matcal_calculus_ode_public_header_tests PRIVATE
        MatCal::Calculus
        MatCal::ODE
    )
    matcal_configure_test_target(matcal_calculus_ode_public_header_tests)
    add_test(NAME matcal_calculus_ode_public_header_tests COMMAND matcal_calculus_ode_public_header_tests)

    add_executable(matcal_calculus_ode_multi_tu_tests
        tests/calculus_ode_multi_tu_main.cpp
        tests/calculus_multi_tu_a.cpp
        tests/ode_multi_tu_b.cpp
    )
    target_link_libraries(matcal_calculus_ode_multi_tu_tests PRIVATE
        MatCal::Calculus
        MatCal::ODE
    )
    matcal_configure_test_target(matcal_calculus_ode_multi_tu_tests)
    add_test(NAME matcal_calculus_ode_multi_tu_tests COMMAND matcal_calculus_ode_multi_tu_tests)

    add_executable(matcal_calculus_only_consumer_tests
        tests/calculus_only_consumer.cpp
    )
    target_link_libraries(matcal_calculus_only_consumer_tests PRIVATE MatCal::Calculus)
    matcal_configure_test_target(matcal_calculus_only_consumer_tests)
    add_test(NAME matcal_calculus_only_consumer_tests COMMAND matcal_calculus_only_consumer_tests)

    add_executable(matcal_ode_only_consumer_tests
        tests/ode_only_consumer.cpp
    )
    target_link_libraries(matcal_ode_only_consumer_tests PRIVATE MatCal::ODE)
    matcal_configure_test_target(matcal_ode_only_consumer_tests)
    add_test(NAME matcal_ode_only_consumer_tests COMMAND matcal_ode_only_consumer_tests)

    matcal_add_nonlinear_test(matcal_nonlinear_tests
        tests/nonlinear_tests.cpp
    )

    matcal_add_least_squares_test(matcal_least_squares_tests
        tests/least_squares_tests.cpp
    )

    add_executable(matcal_nonlinear_least_squares_interpolation_public_header_tests
        tests/nonlinear_least_squares_interpolation_public_headers_main.cpp
        tests/header_nonlinear_self_contained.cpp
        tests/header_least_squares_self_contained.cpp
        tests/header_polynomial_interpolation_self_contained.cpp
    )
    target_link_libraries(matcal_nonlinear_least_squares_interpolation_public_header_tests PRIVATE
        MatCal::Nonlinear
        MatCal::LeastSquares
        MatCal::Interpolation
    )
    matcal_configure_test_target(matcal_nonlinear_least_squares_interpolation_public_header_tests)
    add_test(NAME matcal_nonlinear_least_squares_interpolation_public_header_tests
        COMMAND matcal_nonlinear_least_squares_interpolation_public_header_tests)

    add_executable(matcal_nonlinear_least_squares_interpolation_multi_tu_tests
        tests/nonlinear_least_squares_interpolation_multi_tu_main.cpp
        tests/nonlinear_multi_tu_a.cpp
        tests/least_squares_multi_tu_b.cpp
        tests/interpolation_remaining_multi_tu_c.cpp
    )
    target_link_libraries(matcal_nonlinear_least_squares_interpolation_multi_tu_tests PRIVATE
        MatCal::Nonlinear
        MatCal::LeastSquares
        MatCal::Interpolation
    )
    matcal_configure_test_target(matcal_nonlinear_least_squares_interpolation_multi_tu_tests)
    add_test(NAME matcal_nonlinear_least_squares_interpolation_multi_tu_tests
        COMMAND matcal_nonlinear_least_squares_interpolation_multi_tu_tests)

    add_executable(matcal_nonlinear_only_consumer_tests
        tests/nonlinear_only_consumer.cpp
    )
    target_link_libraries(matcal_nonlinear_only_consumer_tests PRIVATE MatCal::Nonlinear)
    matcal_configure_test_target(matcal_nonlinear_only_consumer_tests)
    add_test(NAME matcal_nonlinear_only_consumer_tests COMMAND matcal_nonlinear_only_consumer_tests)

    add_executable(matcal_least_squares_only_consumer_tests
        tests/least_squares_only_consumer.cpp
    )
    target_link_libraries(matcal_least_squares_only_consumer_tests PRIVATE MatCal::LeastSquares)
    matcal_configure_test_target(matcal_least_squares_only_consumer_tests)
    add_test(NAME matcal_least_squares_only_consumer_tests COMMAND matcal_least_squares_only_consumer_tests)

    add_executable(matcal_interpolation_only_consumer_tests
        tests/interpolation_only_consumer.cpp
    )
    target_link_libraries(matcal_interpolation_only_consumer_tests PRIVATE MatCal::Interpolation)
    matcal_configure_test_target(matcal_interpolation_only_consumer_tests)
    add_test(NAME matcal_interpolation_only_consumer_tests COMMAND matcal_interpolation_only_consumer_tests)

    matcal_add_linalg_test(matcal_linalg_vector_tests
        tests/linalg/vector_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_dense_matrix_tests
        tests/linalg/dense_matrix_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_solver_tests
        tests/linalg/solver_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_scale_contract_tests
        tests/linalg/scale_contract_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_iterative_solver_tests
        tests/linalg/iterative_solver_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_eigen_solver_tests
        tests/linalg/eigen_solver_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_skyline_matrix_tests
        tests/linalg/skyline_matrix_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_skyline_ldlt_tests
        tests/linalg/skyline_ldlt_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_skyline_performance_contract_tests
        tests/linalg/skyline_performance_contract_tests.cpp
    )

    matcal_add_linalg_test(matcal_linalg_public_header_tests
        tests/linalg/public_headers_main.cpp
        tests/linalg/header_vector_self_contained.cpp
        tests/linalg/header_dense_matrix_self_contained.cpp
        tests/linalg/header_solver_types_self_contained.cpp
        tests/linalg/header_dense_solver_self_contained.cpp
        tests/linalg/header_iterative_solvers_self_contained.cpp
        tests/linalg/header_eigen_solvers_self_contained.cpp
        tests/linalg/header_symmetric_skyline_self_contained.cpp
        tests/linalg/header_skyline_ldlt_self_contained.cpp
    )

    matcal_add_linalg_test(matcal_linalg_multi_tu_tests
        tests/linalg/multi_tu_main.cpp
        tests/linalg/multi_tu_a.cpp
        tests/linalg/multi_tu_b.cpp
        tests/linalg/multi_tu_skyline.cpp
        tests/linalg/multi_tu_iterative_eigen.cpp
    )

    add_executable(matcal_multi_tu_link_tests
        tests/multi_tu_main.cpp
        tests/multi_tu_a.cpp
        tests/multi_tu_b.cpp
    )
    target_link_libraries(matcal_multi_tu_link_tests PRIVATE MatCal::Legacy)
    matcal_configure_test_target(matcal_multi_tu_link_tests)
    add_test(NAME matcal_multi_tu_link_tests COMMAND matcal_multi_tu_link_tests)
endif()
