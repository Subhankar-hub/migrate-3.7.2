# Add the source files for the Haru library
set(HARU_SRC_FILES
        libharuc.cc
        libharuc_cnsfonts.cc
        libharuc_cntfonts.cc
        libharuc_jpfonts.cc
        libharuc_krfonts.cc
        libharu_jpeg.cc
        libharu_png.cc
        PdfAnnotations.cc
        PdfAutoPtrObjects.cc
        PdfBorderStyle.cc
        PdfCatalog.cc
        PdfCIDFontDef.cc
        PdfCMap.cc
        PdfCMap_ETen_B5.cc
        PdfCMap_EUC.cc
        PdfCMap_GB_EUC.cc
        PdfCMap_KSC_EUC.cc
        PdfCMap_KSCms_UHC.cc
        PdfCMap_RKSJ.cc
        PdfContents.cc
        PdfDestinations.cc
        PdfDoc.cc
        PdfEncoding.cc
        PdfEncoding_ISO8859.cc
        PdfEncoding_KOI8.cc
        PdfEncoding_MSCP.cc
        PdfEncrypt.cc
        PdfExceptions.cc
        PdfFontDef.cc
        PdfFontDef_Courier.cc
        PdfFontDef_Helvetica.cc
        PdfFontDef_Symbol.cc
        PdfFontDef_Times.cc
        PdfFontDef_ZapfDingbats.cc
        PdfFonts.cc
        PdfInfo.cc
        PdfJpegImage.cc
        PdfList.cc
        PdfLList.cc
        PdfMbFontDef_Batang.cc
        PdfMbFontDef_Dotum.cc
        PdfMbFontDef_Gothic.cc
        PdfMbFontDef_Mincyo.cc
        PdfMbFontDef_Ming.cc
        PdfMbFontDef_SimHei.cc
        PdfMbFontDef_SimSun.cc
        PdfMbFonts.cc
        PdfObjects.cc
        PdfOutlines.cc
        PdfPages.cc
        PdfPngImage.cc
        PdfStreams.cc
        PdfUnicodeText.cc
        PdfUtils.cc
        PdfXObjects.cc
        PdfXref.cc
)

# Add the library
add_library(haru ${HARU_SRC_FILES})