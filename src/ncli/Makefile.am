SRC = nclisr.f90 nclireheat.f90
MOD = nclisr.$(FC_MODEXT) nclireheat.$(FC_MODEXT)

check_PROGRAMS = nclimain
nclimain_SOURCES = $(SRC) nclimain.f90
nclimain_FCFLAGS = -I../$(SRCOMMDIR)
nclimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libncli.la
libncli_la_SOURCES = $(SRC)
libncli_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libncli_la_includedir = $(includedir)/$(SRINCDIR)
libncli_la_include_HEADERS = $(MOD)

man_MANS = aspic_ncli.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
