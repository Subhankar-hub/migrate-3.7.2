/*
 * << H a r u --free pdf library >> -- PdfPngImage.cpp
 *
 * Copyright (c) 1999-2003 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 */

#include "libharu.h"
#include "hpdf_png.h"
#include <png.h>
#include <stdlib.h>
#include <stdio.h>

/*----- PdfPngImage class ----------------------------------------------------*/

void PdfPngImage::LoadFromFile(const char* filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        // Handle error
        return;
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fclose(fp);
        // Handle error
        return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
        fclose(fp);
        // Handle error
        return;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
        fclose(fp);
        // Handle error
        return;
    }

    png_init_io(png_ptr, fp);
    png_read_info(png_ptr, info_ptr);

    if (png_get_bit_depth(png_ptr, info_ptr) == 16) {
        png_set_strip_16(png_ptr);
    }

    if (png_get_color_type(png_ptr, info_ptr) & PNG_COLOR_MASK_ALPHA) {
        png_set_strip_alpha(png_ptr);
    }

    if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png_ptr);
    } else if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_GRAY) {
        png_set_expand_gray_1_2_4_to_8(png_ptr);
    }

    png_bytep* row_pointers = new png_bytep[png_get_image_height(png_ptr, info_ptr)];
    memset(row_pointers, 0x00, png_get_image_height(png_ptr, info_ptr));

    for (int i = 0; i < (int)png_get_image_height(png_ptr, info_ptr); i++) {
        row_pointers[i] = new png_byte[png_get_rowbytes(png_ptr, info_ptr)];
    }

    png_read_image(png_ptr, row_pointers);

    for (int i = 0; i < (int)png_get_image_height(png_ptr, info_ptr); i++) {
        delete[] row_pointers[i];
    }
    delete[] row_pointers;

    fWidth = (unsigned int)png_get_image_width(png_ptr, info_ptr);
    fHeight = (unsigned int)png_get_image_height(png_ptr, info_ptr);
    fBitsPerComponent = (unsigned int)png_get_bit_depth(png_ptr, info_ptr);

    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
    fclose(fp);
}