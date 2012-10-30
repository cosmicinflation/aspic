SRC = ssicommon.f90 ssi1sr.f90 ssi1reheat.f90
MOD = ssi1sr.$(FC_MODEXT) ssi1reheat.$(FC_MODEXT)

check_PROGRAMS = ssi1main
ssi1main_SOURCES = $(SRC) ssi1main.f90
ssi1main_FCFLAGS = -I../$(SRCOMMDIR)
ssi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libssi1.la
libssi1_la_SOURCES = $(SRC)
libssi1_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libssi1_la_includedir = $(includedir)/$(SRINCDIR)
libssi1_la_include_HEADERS = $(MOD)

man_MANS = aspic_ssi1.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



