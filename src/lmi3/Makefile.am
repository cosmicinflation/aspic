SRC = lmi3sr.f90 lmi3reheat.f90
MOD = lmi3sr.$(FC_MODEXT) lmi3reheat.$(FC_MODEXT)

check_PROGRAMS = lmi3main
lmi3main_SOURCES = $(SRC) lmi3main.f90
lmi3main_FCFLAGS = -I../$(SRCOMMDIR)
lmi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liblmi3.la
liblmi3_la_SOURCES = $(SRC)
liblmi3_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblmi3_la_includedir = $(includedir)/$(SRINCDIR)
liblmi3_la_include_HEADERS = $(MOD)

man_MANS = aspic_lmi3.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



