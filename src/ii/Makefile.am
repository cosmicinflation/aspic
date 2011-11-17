SRC = iisr.f90 iireheat.f90
MOD = iisr.$(FC_MODEXT) iireheat.$(FC_MODEXT)

check_PROGRAMS = iimain
iimain_SOURCES = $(SRC) iimain.f90
iimain_FCFLAGS = -I../$(SRCOMMDIR)
iimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libii.la
libii_la_SOURCES = $(SRC)
libii_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libii_la_includedir = $(includedir)/$(SRINCDIR)
libii_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



