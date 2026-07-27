cmake_minimum_required(VERSION 3.20)

project(MatCal
    VERSION 0.1.0
    DESCRIPTION "Legacy C++ numerical calculation library"
    LANGUAGES CXX
)

option(MATCAL_ENABLE_WARNINGS "Enable compiler warnings for MatCal test targets" ON)
option(MATCAL_ENABLE_SANITIZERS "Enable AddressSanitizer/UndefinedBehaviorSanitizer for MatCal test targets when supported" OFF)

add_library(matcal_legacy INTERFACE)
add_library(MatCal::Legacy ALIAS matcal_legacy)

target_compile_features(matcal_legacy INTERFACE cxx_std_20)
target_include_directories(matcal_legacy
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(matcal_legacy PROPERTIES EXPORT_NAME Legacy)

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

include(CTest)

if(BUILD_TESTING)
    add_executable(matcal_characterization_tests
        tests/characterization_tests.cpp
    )
    target_link_libraries(matcal_characterization_tests PRIVATE MatCal::Legacy)
    matcal_configure_test_target(matcal_characterization_tests)
    add_test(NAME matcal_characterization_tests COMMAND matcal_characterization_tests)

    add_executable(matcal_multi_tu_link_tests
        tests/multi_tu_main.cpp
        tests/multi_tu_a.cpp
        tests/multi_tu_b.cpp
    )
    target_link_libraries(matcal_multi_tu_link_tests PRIVATE MatCal::Legacy)
    matcal_configure_test_target(matcal_multi_tu_link_tests)
    add_test(NAME matcal_multi_tu_link_tests COMMAND matcal_multi_tu_link_tests)
endif()
