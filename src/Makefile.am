AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libsrmodels.la
libsrmodels_la_SOURCES =

#forces the linker to use fortran linking
nodist_EXTRA_libsrmodels_la_SOURCES = donotexist.f90

#while the real sources are already built libs
libsrmodels_la_LIBADD =                \
        common/libsrcommon.la          \
        smallf/libsmallf.la            \
        largef/liblargef.la            \
        kksf/libkksf.la                \
        kklt/libkklt.la

#to be build from there
SUBDIRS = common smallf largef kksf kklt

