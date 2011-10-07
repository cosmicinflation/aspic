HYSRC = hysrevol.f90 hyreheat.f90

check_PROGRAMS = hytest
hytest_SOURCES = $(HYSRC) hytest.f90
hytest_FCFLAGS = -I../$(SRCOMMDIR)
hytest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhybrid.la
libhybrid_la_SOURCES = $(HYSRC)
libhybrid_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
