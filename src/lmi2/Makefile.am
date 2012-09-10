SRC = lmi2sr.f90 lmi2reheat.f90
MOD = lmi2sr.$(FC_MODEXT) lmi2reheat.$(FC_MODEXT)

check_PROGRAMS = lmi2main
lmi2main_SOURCES = $(SRC) lmi2main.f90
lmi2main_FCFLAGS = -I../$(SRCOMMDIR)
lmi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblmi2.la
liblmi2_la_SOURCES = $(SRC)
liblmi2_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblmi2_la_includedir = $(includedir)/$(SRINCDIR)
liblmi2_la_include_HEADERS = $(MOD)

man_MANS = aspic_lmi2.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



