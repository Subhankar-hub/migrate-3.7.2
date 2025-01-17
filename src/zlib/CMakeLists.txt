# Add the source files for the zlib library
set(ZLIB_SRC_FILES
        adler32.c
        compress.c
        crc32.c
        deflate.c
        example.c
        gzclose.c
        gzlib.c
        gzread.c
        gzwrite.c
        infback.c
        inffast.c
        inflate.c
        inftrees.c
        minigzip.c
        trees.c
        uncompr.c
        zutil.c
)

# Add the library
add_library(zlib ${ZLIB_SRC_FILES})