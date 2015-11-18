SRC = dicommon.f90 disr.f90

MOD = dicommon.$(FC_MODEXT) 

check_PROGRAMS = dimain
di1main_SOURCES = $(SRC) dimain.f90
di1main_FCFLAGS = -I../$(SRCOMMDIR)
di1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libdi.la
libdi_la_SOURCES = $(SRC)
libdi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libdi_la_includedir = $(includedir)/$(SRINCDIR)
libdi_la_include_HEADERS = $(MOD)

man_MANS = aspic_di.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
