lib_LTLIBRARIES = libsrmodels.la
libsrmodels_la_SOURCES =

#forces the linker to use fortran linking
nodist_EXTRA_libsrmodels_la_SOURCES = donotexist.f90

#while the real sources are already built libs
libsrmodels_la_LIBADD =                \
        common/libsrcommon.la          \
        lfi/liblfi.la                  \
        sfi/libsfi.la                  
#        sfbi/libsfbi.la                \
        bi/libbi.la                    \
        mlfi/libmlfi.la                \
        rmi/librmi.la                  \
        vhi/libvhi.la                 

#to be build from there
SUBDIRS = common lfi sfi
#sfi sfbi bi mlfi rmi vhi

