// PdfCIDWItem.h
#ifndef PDFCIDWITEM_H
#define PDFCIDWITEM_H

class PdfCIDWItem {
public:
    virtual ~PdfCIDWItem();
    // Other member functions and variables
};

class PdfCidWItem1 : public PdfCIDWItem {
public:
    PdfCidWItem1(unsigned int from, unsigned int to, unsigned int value);
    // Other member functions and variables
private:
    unsigned int fFrom;
    unsigned int fTo;
    unsigned int fValue;
};

class PdfCidWItem2 : public PdfCIDWItem {
public:
    PdfCidWItem2(unsigned int from, unsigned int count, unsigned int* values);
    ~PdfCidWItem2();
    // Other member functions and variables
private:
    unsigned int fFrom;
    unsigned int fCount;
    unsigned int* fWidths;
};

#endif // PDFCIDWITEM_H