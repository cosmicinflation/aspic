LARGEFSRC = lfsrevol.f90 lfreheat.f90

check_PROGRAMS = lfrehmain
lfrehmain_SOURCES = $(LARGEFSRC) lfrehmain.f90
lfrehmain_FCFLAGS = -I../$(SRCOMMDIR)
lfrehmain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblargef.la
liblargef_la_SOURCES = $(LARGEFSRC)
liblargef_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)




