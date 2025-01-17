cmake_minimum_required(VERSION 3.10)
project(migrate)

# Set the C and C++ standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# Add subdirectories for different components
add_subdirectory(src)
add_subdirectory(src/haru)
add_subdirectory(src/zlib)

# Include directories
include_directories(src)
include_directories(src/haru)
include_directories(src/zlib)

# Add the main executable
add_executable(migrate src/main.c)
target_link_libraries(migrate zlib haru)