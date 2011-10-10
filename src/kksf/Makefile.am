KKSFSRC = kksfsrevol.f90 kksfreheat.f90

check_PROGRAMS = kksftest
kksftest_SOURCES = $(KKSFSRC) kksftest.f90
kksftest_FCFLAGS = -I../$(SRCOMMDIR) -I../smallf
kksftest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a 
kksftest_LDADD += ../smallf/libsmallf_la-sfreheat.o
kksftest_LDADD += ../smallf/libsmallf_la-sfsrevol.o

noinst_LTLIBRARIES = libkksf.la
libkksf_la_SOURCES = $(KKSFSRC)
libkksf_la_FCFLAGS = -I../$(SRCOMMDIR) -I../smallf $(AM_FCFLAGS)
libkksf_la_includedir = $(includedir)/$(SRINCDIR)
libkksf_la_include_HEADERS = *.$(FC_MODEXT)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
