KKLTSRC = kkltsrevol.f90 kkltreheat.f90

check_PROGRAMS = kklttest
kklttest_SOURCES = $(KKLTSRC) kklttest.f90
kklttest_FCFLAGS = -I../$(SRCOMMDIR)
kklttest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libkklt.la
libkklt_la_SOURCES = $(KKLTSRC)
libkklt_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
