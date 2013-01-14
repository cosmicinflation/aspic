lib_LTLIBRARIES = libaspic.la
libaspic_la_SOURCES =
libaspic_la_LDFLAGS = -version-info @version_info@

libaspic_la_includedir = $(includedir)/$(SRINCDIR)
libaspic_la_include_HEADERS = ../include/aspic.h

#forces the linker to use fortran linking
nodist_EXTRA_libaspic_la_SOURCES = donotexist.f90

#while the real sources are already built libs
libaspic_la_LIBADD =                \
        common/libsrcommon.la          \
        lfi/liblfi.la \
        gmlfi/libgmlfi.la \
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
        lmi/liblmi.la \
        rgi/librgi.la \
	vhi/libvhi.la \
	dsi/libdsi.la \
        ssi/libssi.la \
	oi/liboi.la \
	ai/libai.la \
	cnai/libcnai.la \
	cnbi/libcnbi.la \
	ripi/libripi.la \
	csi/libcsi.la \
	sbi/libsbi.la \
	cnci/libcnci.la \
	rmi/librmi.la \
	lpi/liblpi.la \
	cndi/libcndi.la


#to be build from there
SUBDIRS = common lfi gmlfi rcmi rcqi pni mni esi pli kmii hf1i sfi li ii cwi kmiii hi twi mhi rpi dwi mssmi bsusybi ti bei psni ncki lmi rgi vhi dsi ssi oi ai cnai cnbi ripi csi sbi cnci rmi lpi cndi

