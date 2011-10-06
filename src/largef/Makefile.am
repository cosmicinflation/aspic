LARGEFSRC = lfsrevol.f90 lfreheat.f90

check_PROGRAMS = lftest
lftest_SOURCES = $(LARGEFSRC) lftest.f90
lftest_FCFLAGS = -I../$(SRCOMMDIR)
lftest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblargef.la
liblargef_la_SOURCES = $(LARGEFSRC)
liblargef_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)




