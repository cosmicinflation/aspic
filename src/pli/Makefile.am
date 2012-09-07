SRC = plisr.f90 plireheat.f90
MOD = plisr.$(FC_MODEXT) plireheat.$(FC_MODEXT)

check_PROGRAMS = plimain
plimain_SOURCES = $(SRC) plimain.f90
plimain_FCFLAGS = -I../$(SRCOMMDIR)
plimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libpli.la
libpli_la_SOURCES = $(SRC)
libpli_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libpli_la_includedir = $(includedir)/$(SRINCDIR)
libpli_la_include_HEADERS = $(MOD)

man_MANS = libpli.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



