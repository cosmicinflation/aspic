SRC = fisr.f90 fireheat.f90
MOD = fisr.$(FC_MODEXT) fireheat.$(FC_MODEXT)

check_PROGRAMS = fimain
fimain_SOURCES = $(SRC) fimain.f90
fimain_FCFLAGS = -I../$(SRCOMMDIR)
fimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libfi.la
libfi_la_SOURCES = $(SRC)
libfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libfi_la_includedir = $(includedir)/$(SRINCDIR)
libfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_fi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
