lib_LTLIBRARIES = libaspic@PREC_SUFFIX@.la
libaspic@PREC_SUFFIX@_la_SOURCES =
libaspic@PREC_SUFFIX@_la_LDFLAGS = -version-info @version_info@

libaspic@PREC_SUFFIX@_la_includedir = $(includedir)/$(SRINCDIR)
libaspic@PREC_SUFFIX@_la_include_HEADERS = ../include/aspic.h

#forces the linker to use fortran linking
nodist_EXTRA_libaspic@PREC_SUFFIX@_la_SOURCES = donotexist.f90


#while the real sources are already built libs
libaspic@PREC_SUFFIX@_la_LIBADD =                \
        common/libsrcommon.la          \
        lfi/liblfi.la \
        gmlfi/libgmlfi.la \
        rcmi/librcmi.la \
	rcqi/librcqi.la \
	ni/libni.la \
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
        gmssmi/libgmssmi.la \
        bsusybi/libbsusybi.la \
        ti/libti.la \
        bei/libbei.la \
        psni/libpsni.la \
        ncki/libncki.la \
        lmi/liblmi.la \
        rgi/librgi.la \
	vhi/libvhi.la \
	dsi/libdsi.la \
        ssbi/libssbi.la \
	oi/liboi.la \
	ai/libai.la \
	cnai/libcnai.la \
	cnbi/libcnbi.la \
	oripi/liboripi.la \
	csi/libcsi.la \
	sbi/libsbi.la \
	cnci/libcnci.la \
	rmi/librmi.la \
	lpi/liblpi.la \
	cndi/libcndi.la \
	rchi/librchi.la \
        osti/libosti.la \
        wri/libwri.la \
        imi/libimi.la \
        gripi/libgripi.la \
        bi/libbi.la \
        kklti/libkklti

#to be build from there
SUBDIRS = common lfi gmlfi rcmi rcqi ni esi pli kmii hf1i sfi li ii cwi kmiii hi twi \
          mhi rpi dwi gmssmi bsusybi ti bei psni ncki lmi rgi vhi dsi ssbi oi ai cnai cnbi \
          oripi csi sbi cnci rmi lpi cndi rchi osti wri imi gripi bi kklti

#make clean does not remove modules by default
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
