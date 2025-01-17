# Add the source files for the main project
set(SRC_FILES
        aic.c
        altivec.c
        autotune.c
        bayes.c
        broyden.c
        calculator.c
        combroyden2.c
        correlation.c
        data.c
        definitions.h
        fst.c
        gammalike.c
        heating.c
        histogram.c
        inthp.c
        joint-chains.c
        laguerre.c
        lrt.c
        main.c
        mcmc1.c
        mcmc2.c
        mcmc2.h
        menu.c
        migevents.c
        migrate_mpi.c
        mutationmodel.c
        options.c
        pretty.c
        priors.c
        profile.c
        random.c
        reporter.c
        sequence.c
        sighandler.c
        skyline.c
        slice.c
        sort.c
        spline.c
        tools.c
        tree.c
        uep.c
        watterson.c
        world.c
        znzlib.c
)

# Add the library
add_library(migrate_lib ${SRC_FILES})

# Link with zlib and haru libraries
target_link_libraries(migrate_lib zlib haru)