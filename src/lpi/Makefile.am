SRC = lpisr.f90 lpireheat.f90
MOD = lpisr.$(FC_MODEXT) lpireheat.$(FC_MODEXT)

check_PROGRAMS = lpimain
lpimain_SOURCES = $(SRC) lpimain.f90
lpimain_FCFLAGS = -I../$(SRCOMMDIR)
lpimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
mlfimain_SOURCES = $(SRC) mlfimain.f90
mlfimain_FCFLAGS = -I../$(SRCOMMDIR)
mlfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblpi.la
liblpi_la_SOURCES = $(SRC)
liblpi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblpi_la_includedir = $(includedir)/$(SRINCDIR)
liblpi_la_include_HEADERS = $(MOD)

man_MANS = aspic_lpi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



