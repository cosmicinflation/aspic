SRC = rcipicommon.f90 rcipisr.f90 rcipireheat.f90
MOD = rcipicommon.$(FC_MODEXT) rcipisr.$(FC_MODEXT) rcipireheat.$(FC_MODEXT)

check_PROGRAMS = rcipimain
rcipimain_SOURCES = $(SRC) rcipimain.f90
rcipimain_FCFLAGS = -I../$(SRCOMMDIR)
rcipimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librcipi.la
librcipi_la_SOURCES = $(SRC)
librcipi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpi_la_includedir = $(includedir)/$(SRINCDIR)
librpi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rcipi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
