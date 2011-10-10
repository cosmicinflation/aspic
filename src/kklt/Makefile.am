KKLTSRC = kkltsrevol.f90 kkltreheat.f90
KKLTMOD =  kkltsrevol.$(FC_MODEXT) kkltreheat.$(FC_MODEXT)

check_PROGRAMS = kklttest
kklttest_SOURCES = $(KKLTSRC) kklttest.f90
kklttest_FCFLAGS = -I../$(SRCOMMDIR)
kklttest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libkklt.la
libkklt_la_SOURCES = $(KKLTSRC)
libkklt_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libkklt_la_includedir = $(includedir)/$(SRINCDIR)
libkklt_la_include_HEADERS = $(KKLTMOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
