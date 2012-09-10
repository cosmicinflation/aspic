SRC = tisr.f90 tireheat.f90
MOD = tisr.$(FC_MODEXT) tireheat.$(FC_MODEXT)

check_PROGRAMS = timain
timain_SOURCES = $(SRC) timain.f90
timain_FCFLAGS = -I../$(SRCOMMDIR)
timain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libti.la
libti_la_SOURCES = $(SRC)
libti_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libti_la_includedir = $(includedir)/$(SRINCDIR)
libti_la_include_HEADERS = $(MOD)

man_MANS = aspic_ti.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
