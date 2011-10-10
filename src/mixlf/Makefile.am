MIXLFSRC = mixlfsrevol.f90 mixlfreheat.f90
MIXLFMOD = mixlfsrevol.$(FC_MODEXT) mixlfreheat.$(FC_MODEXT)

check_PROGRAMS = mixlftest
mixlftest_SOURCES = $(MIXLFSRC) mixlftest.f90
mixlftest_FCFLAGS = -I../$(SRCOMMDIR)
mixlftest_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libmixlf.la
libmixlf_la_SOURCES = $(MIXLFSRC)
libmixlf_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libmixlf_la_includedir = $(includedir)/$(SRINCDIR)
libmixlf_la_include_HEADERS = $(MIXLFMOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
