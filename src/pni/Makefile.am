SRC = pnisr.f90 pnireheat.f90
MOD = pnisr.$(FC_MODEXT) pnireheat.$(FC_MODEXT)

check_PROGRAMS = pnimain
pnimain_SOURCES = $(SRC) pnimain.f90
pnimain_FCFLAGS = -I../$(SRCOMMDIR)
pnimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libpni.la
libpni_la_SOURCES = $(SRC)
libpni_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libpni_la_includedir = $(includedir)/$(SRINCDIR)
libpni_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



