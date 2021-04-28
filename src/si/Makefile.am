SRC = sisr.f90 sireheat.f90
MOD = sisr.$(FC_MODEXT) sireheat.$(FC_MODEXT)

check_PROGRAMS = simain
simain_SOURCES = $(SRC) simain.f90
simain_FCFLAGS = -I../$(SRCOMMDIR)
simain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsi.la
libsi_la_SOURCES = $(SRC)
libsi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsi_la_includedir = $(includedir)/$(SRINCDIR)
libsi_la_include_HEADERS = $(MOD)

man_MANS = aspic_si.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
