SRC = displine.f90 dicommon.f90 disr.f90

MOD = displine.$(FC_MODEXT) dicommon.$(FC_MODEXT) disr.$(FC_MODEXT)

check_PROGRAMS = dimain
dimain_SOURCES = $(SRC) dimain.f90
dimain_FCFLAGS = -I../$(SRCOMMDIR)
dimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

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
