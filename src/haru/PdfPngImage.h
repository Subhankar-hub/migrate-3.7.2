#ifndef PDFPNGIMAGE_H
#define PDFPNGIMAGE_H

#include <png.h>

class PdfPngImage {
public:
    PdfPngImage();
    ~PdfPngImage();
    void LoadFromFile(const char* filename);

private:
    unsigned int fWidth;
    unsigned int fHeight;
    unsigned int fBitsPerComponent;
};

#endif // PDFPNGIMAGE_H