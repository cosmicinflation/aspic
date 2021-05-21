SRC = hispline.F90 hicommon.f90 hisr.f90 hireheat.f90

MOD = hispline.$(FC_MODEXT) hicommon.$(FC_MODEXT) hisr.$(FC_MODEXT) hireheat.$(FC_MODEXT)

check_PROGRAMS = himain
himain_SOURCES = $(SRC) himain.f90
himain_FCFLAGS = -I../$(SRCOMMHIR)
himain_LDADD = ../$(SRCOMMHIR)/libsrcommon.a

noinst_LTLIBRARIES = libhi.la
libhi_la_SOURCES = $(SRC)
libhi_la_FCFLAGS = -I../$(SRCOMMHIR) $(AM_FCFLAGS)
libhi_la_includehir = $(includehir)/$(SRINCHIR)
libhi_la_include_HEADERS = $(MOD)

man_MANS = aspic_hi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
