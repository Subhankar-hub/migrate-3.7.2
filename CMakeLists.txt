cmake_minimum_required(VERSION 3.10)
project(migrate)

# Set the C and C++ standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES
        src/aic.c
        src/altivec.c
        src/autotune.c
        src/bayes.c
        src/broyden.c
        src/calculator.c
        src/combroyden2.c
        src/correlation.c
        src/data.c
        src/definitions.h
        src/fst.c
        src/gammalike.c
        src/heating.c
        src/histogram.c
        src/inthp.c
        src/joint-chains.c
        src/laguerre.c
        src/lrt.c
        src/main.c
        src/mcmc1.c
        src/mcmc2.c
        src/mcmc2.h
        src/menu.c
        src/migevents.c
        src/migrate_mpi.c
        src/mutationmodel.c
        src/options.c
        src/pretty.c
        src/priors.c
        src/profile.c
        src/random.c
        src/reporter.c
        src/sequence.c
        src/sighandler.c
        src/skyline.c
        src/slice.c
        src/sort.c
        src/spline.c
        src/tools.c
        src/tree.c
        src/uep.c
        src/watterson.c
        src/world.c
        src/znzlib.c
)
# Add subdirectories for different components
add_subdirectory(src)
add_subdirectory(src/haru)
add_subdirectory(src/zlib)

# Include directories
include_directories(src)
include_directories(src/haru)
include_directories(src/zlib)
add_library(migrate_lib STATIC ${SOURCE_FILES})
# Add the main executable
add_executable(migrate src/main.c)
target_link_libraries(migrate_lib -lz -lpng -ljpeg)
target_link_libraries(migrate migrate_lib zlib haru)