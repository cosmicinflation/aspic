lib_LTLIBRARIES = libaspic@PREC_SUFFIX@.la
libaspic@PREC_SUFFIX@_la_SOURCES =
libaspic@PREC_SUFFIX@_la_LDFLAGS = -version-info @version_info@

libaspic@PREC_SUFFIX@_la_includedir = $(includedir)/$(SRINCDIR)
libaspic@PREC_SUFFIX@_la_include_HEADERS = ../include/aspic@PREC_SUFFIX@.h

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
        kklti/libkklti.la \
        nfi/libnfi.la \
	ccsi/libccsi.la \
	di/libdi.la \
	hni/libhni.la \
	ncli/libncli.la \
	vfmi/libvfmi.la \
        fi/libfi.la \
	ahi/libahi.la \
	sbki/libsbki.la \
	sdi/libsdi.la \
	saai/libsaai.la \
	sabi/libsabi.la \
	saci/libsaci.la \
	hbi/libhbi.la \
	shi/libshi.la \
	rpqdi/librpqdi.la \
	rpqti/librpqti.la \
	rcpi/librcpi.la \
	gdwi/libgdwi.la \
	saii/libsaii.la \
	saiii/libsaiii.la	

#to be build from there
SUBDIRS = lfi gmlfi rcmi rcqi ni esi pli kmii hf1i sfi li ii cwi kmiii hi twi \
          mhi rpi dwi gmssmi bsusybi ti bei psni ncki lmi rgi vhi dsi ssbi oi \
          ai cnai cnbi oripi csi sbi cnci rmi lpi cndi rchi osti wri imi gripi \
          bi kklti nfi ccsi di hni ncli vfmi fi ahi sbki sdi saai sabi saci hbi \
          shi rpqdi rpqti rcpi gdwi saii saiii

#in parallel (bugged)
#.PHONY: #(SUBDIRS)
##(SUBDIRS):
#	#(MAKE) -C #@ #(MAKECMDGOALS)
#all install check test: #(SUBDIRS)

#make clean does not remove modules by default
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)


.NOTPARALLEL:
