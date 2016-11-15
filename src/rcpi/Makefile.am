SRC = rcpicommon.f90 rcpisr.f90 rcpireheat.f90
MOD = rcpicommon.$(FC_MODEXT) rcpisr.$(FC_MODEXT) rcpireheat.$(FC_MODEXT)

check_PROGRAMS = rcpimain
rcpimain_SOURCES = $(SRC) rcpimain.f90
rcpimain_FCFLAGS = -I../$(SRCOMMDIR)
rcpimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librcpi.la
librcpi_la_SOURCES = $(SRC)
librcpi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpi_la_includedir = $(includedir)/$(SRINCDIR)
librpi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rcpi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
