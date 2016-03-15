SRC = abieos.f90 abisr.f90 abireheat.f90

MOD = abieos.$(FC_MODEXT) abisr.$(FC_MODEXT) abireheat.$(FC_MODEXT)

check_PROGRAMS = abimain
abimain_SOURCES = $(SRC) abimain.f90
abimain_FCFLAGS = -I../$(SRCOMMDIR)
abimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libabi.la
libabi_la_SOURCES = $(SRC)
libabi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libabi_la_includedir = $(includedir)/$(SRINCDIR)
libabi_la_include_HEADERS = $(MOD)

man_MANS = aspic_abi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
