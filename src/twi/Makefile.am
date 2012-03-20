SRC = twisr.f90 twireheat.f90
MOD = twisr.$(FC_MODEXT) twireheat.$(FC_MODEXT)

check_PROGRAMS = twimain
twimain_SOURCES = $(SRC) twimain.f90
twimain_FCFLAGS = -I../$(SRCOMMDIR)
twimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libtwi.la
libtwi_la_SOURCES = $(SRC)
libtwi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libtwi_la_includedir = $(includedir)/$(SRINCDIR)
libtwi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



