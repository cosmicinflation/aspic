SRC = rpicommon.f90 rpi1sr.f90 rpi1reheat.f90 rpi2sr.f90 rpi2reheat.f90
MOD = rpicommon.$(FC_MODEXT) rpi1sr.$(FC_MODEXT) rpi1reheat.$(FC_MODEXT) \
	rpi2sr.$(FC_MODEXT) rpi2reheat.$(FC_MODEXT)

check_PROGRAMS = rpi1main rpi2main
rpi1main_SOURCES = $(SRC) rpi1main.f90
rpi1main_FCFLAGS = -I../$(SRCOMMDIR)
rpi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

rpi2main_SOURCES = $(SRC) rpi2main.f90
rpi2main_FCFLAGS = -I../$(SRCOMMDIR)
rpi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librpi.la
librpi_la_SOURCES = $(SRC)
librpi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpi_la_includedir = $(includedir)/$(SRINCDIR)
librpi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rpi1.3 aspic_rpi2.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



