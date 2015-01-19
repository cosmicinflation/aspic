SRC = ccsicommon.f90 ccsi1sr.f90 ccsi1reheat.f90

MOD = ccsicommon.$(FC_MODEXT) ccsi1sr.$(FC_MODEXT) ccsi1reheat.$(FC_MODEXT)

check_PROGRAMS = ccsi1main
ccsi1main_SOURCES = $(SRC) ccsi1main.f90
ccsi1main_FCFLAGS = -I../$(SRCOMMDIR)
ccsi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a


noinst_LTLIBRARIES = libccsi.la
libccsi_la_SOURCES = $(SRC)
libccsi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libccsi_la_includedir = $(includedir)/$(SRINCDIR)
libccsi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ccsi1.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
