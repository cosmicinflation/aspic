SRC = gdwisr.f90 gdwireheat.f90
MOD = gdwisr.$(FC_MODEXT) gdwireheat.$(FC_MODEXT)

check_PROGRAMS = gdwimain
gdwimain_SOURCES = $(SRC) gdwimain.f90
gdwimain_FCFLAGS = -I../$(SRCOMMDIR)
gdwimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libgdwi.la
libgdwi_la_SOURCES = $(SRC)
libgdwi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libgdwi_la_includedir = $(includedir)/$(SRINCDIR)
libgdwi_la_include_HEADERS = $(MOD)

man_MANS = aspic_gdwi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
