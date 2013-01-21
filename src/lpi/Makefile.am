SRC = lpicommon.f90 lpi1sr.f90 lpi1reheat.f90 lpi2sr.f90 \
	lpi2reheat.f90 lpi3sr.f90 lpi3reheat.f90
MOD = lpicommon.$(FC_MODEXT) lpi1sr.$(FC_MODEXT) lpi1reheat.$(FC_MODEXT) \
	lpi2sr.$(FC_MODEXT) lpi2reheat.$(FC_MODEXT) lpi3sr.$(FC_MODEXT) \
	lpi3reheat.$(FC_MODEXT)

check_PROGRAMS = lpi1main lpi2main lpi3main
lpi1main_SOURCES = $(SRC) lpi1main.f90
lpi1main_FCFLAGS = -I../$(SRCOMMDIR)
lpi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
lpi2main_SOURCES = $(SRC) lpi2main.f90
lpi2main_FCFLAGS = -I../$(SRCOMMDIR)
lpi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
lpi3main_SOURCES = $(SRC) lpi3main.f90
lpi3main_FCFLAGS = -I../$(SRCOMMDIR)
lpi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a


noinst_LTLIBRARIES = liblpi.la
liblpi_la_SOURCES = $(SRC)
liblpi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liblpi_la_includedir = $(includedir)/$(SRINCDIR)
liblpi_la_include_HEADERS = $(MOD)

man_MANS = aspic_lpi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
