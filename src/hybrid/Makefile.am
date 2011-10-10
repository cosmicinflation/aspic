HYSRC = hysrevol.f90 hyreheat.f90
HYMOD = hysrevol.$(FC_MODEXT) hyreheat.$(FC_MODEXT)

check_PROGRAMS = hytest
hytest_SOURCES = $(HYSRC) hytest.f90
hytest_FCFLAGS = -I../$(SRCOMMDIR)
hytest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhybrid.la
libhybrid_la_SOURCES = $(HYSRC)
libhybrid_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhybrid_la_includedir = $(includedir)/$(SRINCDIR)
libhybrid_la_include_HEADERS = $(HYMOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
