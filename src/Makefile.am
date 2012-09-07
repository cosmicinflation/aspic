lib_LTLIBRARIES = libaspic.la
libaspic_la_SOURCES =
libascpi_la_LDFLAGS = -version-info @version_info@

#forces the linker to use fortran linking
nodist_EXTRA_libsrmodels_la_SOURCES = donotexist.f90

#while the real sources are already built libs
libaspic_la_LIBADD =                \
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
	kmiii/libkmiii.la \
	hi/libhi.la \
        twi/libtwi.la \
        mhi/libmhi.la \
        rpi/librpi.la \
        dwi/libdwi.la \
        mssmi/libmssmi.la \
        bsusybi/libbsusybi.la \
        ti/libti.la \
        bei/libbei.la \
        psni/libpsni.la \
        ncki/libncki.la \
        lmi1/liblmi1.la \
        lmi2/liblmi2.la \
        lmi3/liblmi3.la
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
SUBDIRS = common lfi mlfi rcmi rcqi pni mni esi pli kmii hf1i sfi li ii cwi kmiii hi twi mhi rpi dwi mssmi bsusybi ti bei psni ncki lmi1 lmi2 lmi3

