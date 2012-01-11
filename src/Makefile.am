lib_LTLIBRARIES = libsrmodels.la
libsrmodels_la_SOURCES =
libsrmodels_la_LDFLAGS = -version-info @version_info@

#forces the linker to use fortran linking
nodist_EXTRA_libsrmodels_la_SOURCES = donotexist.f90

#while the real sources are already built libs
libsrmodels_la_LIBADD =                \
        common/libsrcommon.la          \
        lfi/liblfi.la \
        mlfi/libmlfi.la \
        rcmi/librcmi.la \
	rcqi/librcqi.la \
	pni/libpni.la \
	mni/libmni.la \
	esi/libesi.la \
	pli/libpli.la \
	kmii/libkmii.la \
	hf1i/libhf1i.la \
	sfi/libsfi.la \
	li/libli.la \
	ii/libii.la \
	cwi/libcwi.la \
	kmiii/libkmiii.la
#	bi/libbi.la \
	sfbi/libsfbi.la \
	kmiii/libkmiii.la \
	cwi/libcwi.la \
	hf2i/libhf2i.la \
	rmi/librmi.la \
	vhi/libvhi.la \
	dsi/libdsi.la \
	gmi/libgmi.la

#to be build from there
SUBDIRS = common lfi mlfi rcmi rcqi pni mni esi pli kmii hf1i sfi li ii cwi kmiii
# sfbi bi hf2i rmi vhi dsi gmi

