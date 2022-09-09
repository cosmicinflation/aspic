SRC = ccsicommon.F90 ccsicomreh.f90 ccsi1sr.f90 ccsi1reheat.f90 ccsi2sr.f90 ccsi2reheat.f90 \
ccsi3sr.F90 ccsi3reheat.f90

MOD = ccsicommon.$(FC_MODEXT) ccsicomreh.$(FC_MODEXT) \
	ccsi1sr.$(FC_MODEXT) ccsi1reheat.$(FC_MODEXT) \
	ccsi2sr.$(FC_MODEXT) ccsi2reheat.$(FC_MODEXT) \
	ccsi3sr.$(FC_MODEXT) ccsi3reheat.$(FC_MODEXT)

check_PROGRAMS = ccsi1main ccsi2main ccsi3main
ccsi1main_SOURCES = $(SRC) ccsi1main.f90
ccsi1main_FCFLAGS = -I../$(SRCOMMDIR)
ccsi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

ccsi2main_SOURCES = $(SRC) ccsi2main.f90
ccsi2main_FCFLAGS = -I../$(SRCOMMDIR)
ccsi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

ccsi3main_SOURCES = $(SRC) ccsi3main.f90
ccsi3main_FCFLAGS = -I../$(SRCOMMDIR)
ccsi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libccsi.la
libccsi_la_SOURCES = $(SRC)
libccsi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libccsi_la_includedir = $(includedir)/$(SRINCDIR)
libccsi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ccsi1.3 aspic_ccsi2.3 aspic_ccsi3.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
