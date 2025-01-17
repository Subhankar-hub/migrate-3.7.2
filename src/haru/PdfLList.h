#ifndef PDFLLIST_H
#define PDFLLIST_H

class PdfLList; // Forward declaration

class PdfLNode {
public:
    PdfLNode();
    virtual ~PdfLNode();

    PdfLNode* GetNext() const { return fNext; }
    PdfLNode* GetPrev() const { return fPrev; }

private:
    PdfLNode* fPrev;
    PdfLNode* fNext;
    PdfLList* fParent;
    
    friend class PdfLList;
};

class PdfLList {
public:
    PdfLList();
    virtual ~PdfLList();

    bool Add(PdfLNode* node, PdfLNode* atNode);
    bool Unlink(PdfLNode* node);
    int CountItems();
    void Clear();
    bool IsEmpty() const { return fFirstNode == nullptr; }

private:
    PdfLNode* fFirstNode;
    PdfLNode* fLastNode;
};

#endif // PDFLLIST_H