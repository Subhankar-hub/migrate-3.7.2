#ifndef HPDF_PNG_H
#define HPDF_PNG_H

#include <png.h>
#include "libharu.h"
#include "hpdf_types.h" // Ensure this header is included for HPDF_STATUS and related types

#ifdef __cplusplus
extern "C" {
#endif

// Function to load PNG image from file
HPDF_STATUS HPDF_LoadPngImageFromFile(PdfPngImage *image, const char *filename);

// Function to load PNG image from memory
HPDF_STATUS HPDF_LoadPngImageFromMemory(PdfPngImage *image, const HPDF_BYTE *buffer, HPDF_UINT size);

#ifdef __cplusplus
}
#endif

#endif // HPDF_PNG_H