SMALLFSRC = sfsrevol.f90 sfreheat.f90

check_PROGRAMS = sftest
sftest_SOURCES = $(SMALLFSRC) sftest.f90
sftest_FCFLAGS = -I../$(SRCOMMDIR)
sftest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsmallf.la
libsmallf_la_SOURCES = $(SMALLFSRC)
libsmallf_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
