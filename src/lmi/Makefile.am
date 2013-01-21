SRC = lmicommon.f90 lmi1sr.f90 lmi1reheat.f90 lmi2sr.f90 lmi2reheat.f90
MOD = lmicommon.$(FC_MODEXT) lmi1sr.$(FC_MODEXT) lmi1reheat.$(FC_MODEXT) \
	lmi2sr.$(FC_MODEXT) lmi2reheat.$(FC_MODEXT)

check_PROGRAMS = lmi1main lmi2main
lmi1main_SOURCES = $(SRC) lmi1main.f90
lmi1main_FCFLAGS = -I../$(SRCOMMDIR)
lmi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
lmi2main_SOURCES = $(SRC) lmi2main.f90
lmi2main_FCFLAGS = -I../$(SRCOMMDIR)
lmi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblmi.la
liblmi_la_SOURCES = $(SRC)
liblmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblmi_la_includedir = $(includedir)/$(SRINCDIR)
liblmi_la_include_HEADERS = $(MOD)

man_MANS = aspic_lmi1.3 aspic_lmi2.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
