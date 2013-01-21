SRC = dwisr.f90 dwireheat.f90
MOD = dwisr.$(FC_MODEXT) dwireheat.$(FC_MODEXT)

check_PROGRAMS = dwimain
dwimain_SOURCES = $(SRC) dwimain.f90
dwimain_FCFLAGS = -I../$(SRCOMMDIR)
dwimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libdwi.la
libdwi_la_SOURCES = $(SRC)
libdwi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libdwi_la_includedir = $(includedir)/$(SRINCDIR)
libdwi_la_include_HEADERS = $(MOD)

man_MANS = aspic_dwi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
