SRC = vhisr.f90 vhireheat.f90
MOD = vhisr.$(FC_MODEXT) vhireheat.$(FC_MODEXT)

check_PROGRAMS = vhimain
vhimain_SOURCES = $(SRC) vhimain.f90
vhimain_FCFLAGS = -I../$(SRCOMMDIR)
vhimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libvhi.la
libvhi_la_SOURCES = $(SRC)
libvhi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libvhi_la_includedir = $(includedir)/$(SRINCDIR)
libvhi_la_include_HEADERS = $(MOD)

man_MANS = aspic_vhi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
