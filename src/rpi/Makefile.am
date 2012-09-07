SRC = rpisr.f90 rpireheat.f90
MOD = rpisr.$(FC_MODEXT) rpireheat.$(FC_MODEXT)

check_PROGRAMS = rpimain
rpimain_SOURCES = $(SRC) rpimain.f90
rpimain_FCFLAGS = -I../$(SRCOMMDIR)
rpimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librpi.la
librpi_la_SOURCES = $(SRC)
librpi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpi_la_includedir = $(includedir)/$(SRINCDIR)
librpi_la_include_HEADERS = $(MOD)

man_MANS = librpi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



