INCLUDEPATH += $$system(INCDIR)
INCLUDEPATH +=  /opt/root/install/include

HEADERS += \
    StRoot/StV0Maker/StV0Type.h \
    StRoot/StV0Maker/StV0Maker.h \
    StRoot/StV0Maker/StV0Dst.h \
    StRoot/StV0Maker/StDcaService.h

SOURCES += \
    temp_gccflags.c \
    StRoot/StV0Maker/StV0Maker.cxx \
    StRoot/StV0Maker/StDcaService.cxx \
    reco.C

DISTFILES += \
    input/Scheduler.xml
