cmake_minimum_required(VERSION 3.20)
project(caliper VERSION 1.4.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/W4 /WX /O2)
else()
    add_compile_options(-Wall -Wextra -Werror -O3 -pedantic)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/c)

# Caliper Native Shared Library
add_library(caliper SHARED
    src/c/caliper_io.c
    src/c/caliper_bridge.c
)
target_link_libraries(caliper PRIVATE m)

# Caliper CLI Executable
add_executable(caliper-cli
    src/c/caliper_io.c
    src/c/caliper_bridge.c
    src/c/caliper_main.c
)
target_link_libraries(caliper-cli PRIVATE m)

# C Tests
enable_testing()

add_executable(test_c_io
    src/c/caliper_io.c
    tests/c/test_c_io.c
)
target_link_libraries(test_c_io PRIVATE m)
add_test(NAME TestCaliperIO COMMAND test_c_io)

add_executable(test_c_bridge
    src/c/caliper_io.c
    src/c/caliper_bridge.c
    tests/c/test_c_bridge.c
)
target_link_libraries(test_c_bridge PRIVATE m)
add_test(NAME TestCaliperBridge COMMAND test_c_bridge)
