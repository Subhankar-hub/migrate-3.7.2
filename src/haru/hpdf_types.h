#ifndef HPDF_TYPES_H
#define HPDF_TYPES_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Basic types
typedef uint8_t HPDF_BYTE;
typedef int32_t HPDF_INT;
typedef uint32_t HPDF_UINT;
typedef uint32_t HPDF_STATUS;

// Boolean type
typedef int HPDF_BOOL;
#define HPDF_TRUE 1
#define HPDF_FALSE 0

// Error codes
#define HPDF_OK 0
#define HPDF_ERROR 1

// Declare PdfPngImage as a placeholder, since it is defined elsewhere
class PdfPngImage;

#ifdef __cplusplus
}
#endif

#endif // HPDF_TYPES_H