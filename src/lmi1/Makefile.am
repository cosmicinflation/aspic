SRC = lmi1sr.f90 lmi1reheat.f90
MOD = lmi1sr.$(FC_MODEXT) lmi1reheat.$(FC_MODEXT)

check_PROGRAMS = lmi1main
lmi1main_SOURCES = $(SRC) lmi1main.f90
lmi1main_FCFLAGS = -I../$(SRCOMMDIR)
lmi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblmi1.la
liblmi1_la_SOURCES = $(SRC)
liblmi1_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblmi1_la_includedir = $(includedir)/$(SRINCDIR)
liblmi1_la_include_HEADERS = $(MOD)

man_MANS = liblmi1.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



