SRC = nmlficommon.f90 nmlficomreh.f90 nmlfi1sr.f90 nmlfi1reheat.f90 \
nmlfi2sr.f90 nmlfi2reheat.f90 nmlfi3sr.f90 nmlfi3reheat.f90

MOD = nmlficommon.$(FC_MODEXT) nmlficomreh.$(FC_MODEXT) \
	nmlfi1sr.$(FC_MODEXT) nmlfi1reheat.$(FC_MODEXT) \
	nmlfi2sr.$(FC_MODEXT) nmlfi2reheat.$(FC_MODEXT) \
	nmlfi3sr.$(FC_MODEXT) nmlfi3reheat.$(FC_MODEXT)

check_PROGRAMS = nmlfi1main nmlfi2main nmlfi3main
nmlfi1main_SOURCES = $(SRC) nmlfi1main.f90
nmlfi1main_FCFLAGS = -I../$(SRCOMMDIR)
nmlfi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

nmlfi2main_SOURCES = $(SRC) nmlfi2main.f90
nmlfi2main_FCFLAGS = -I../$(SRCOMMDIR)
nmlfi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

nmlfi3main_SOURCES = $(SRC) nmlfi3main.f90
nmlfi3main_FCFLAGS = -I../$(SRCOMMDIR)
nmlfi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libnmlfi.la
libnmlfi_la_SOURCES = $(SRC)
libnmlfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libnmlfi_la_includedir = $(includedir)/$(SRINCDIR)
libnmlfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_nmlfi1.3 aspic_nmlfi2.3 aspic_nmlfi3.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
