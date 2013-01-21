SRC = lfisr.f90 lfireheat.f90
MOD = lfisr.$(FC_MODEXT) lfireheat.$(FC_MODEXT)

check_PROGRAMS = lfimain
lfimain_SOURCES = $(SRC) lfimain.f90
lfimain_FCFLAGS = -I../$(SRCOMMDIR)
lfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblfi.la
liblfi_la_SOURCES = $(SRC)
liblfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblfi_la_includedir = $(includedir)/$(SRINCDIR)
liblfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_lfi.3


clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)

.NOTPARALLEL:
