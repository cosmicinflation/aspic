SRC = nmlficommon.f90 nmlfisr.f90 nmlfireheat.f90

MOD = nmlficommon.$(FC_MODEXT) nmlfisr.$(FC_MODEXT) nmlfireheat.$(FC_MODEXT)

check_PROGRAMS = nmlfimain
nmlfimain_SOURCES = $(SRC) nmlfimain.f90
nmlfimain_FCFLAGS = -I../$(SRCOMMDIR)
nmlfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libnmlfi.la
libnmlfi_la_SOURCES = $(SRC)
libnmlfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libnmlfi_la_includedir = $(includedir)/$(SRINCDIR)
libnmlfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_nmlfi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
