RMSRC = rmsrevol.f90 rmreheat.f90

check_PROGRAMS = rmtest
rmtest_SOURCES = $(RMSRC) rmtest.f90
rmtest_FCFLAGS = -I../$(SRCOMMDIR)
rmtest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librmass.la
librmass_la_SOURCES = $(RMSRC)
librmass_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librmass_la_includedir = $(includedir)/$(SRINCDIR)
librmass_la_include_HEADERS = *.$(FC_MODEXT)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
