SRC = esisr.f90 esireheat.f90
MOD = esisr.$(FC_MODEXT) esireheat.$(FC_MODEXT)

check_PROGRAMS = esimain
esimain_SOURCES = $(SRC) esimain.f90
esimain_FCFLAGS = -I../$(SRCOMMDIR)
esimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libesi.la
libesi_la_SOURCES = $(SRC)
libesi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libesi_la_includedir = $(includedir)/$(SRINCDIR)
libesi_la_include_HEADERS = $(MOD)

man_MANS = aspic_esi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



