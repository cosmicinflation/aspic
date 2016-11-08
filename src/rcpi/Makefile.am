SRC = rcpicommon.f90 rcp1isr.f90 
MOD = rcpicommon.$(FC_MODEXT) rcp1isr.$(FC_MODEXT)

check_PROGRAMS = rcpi1main
rpqdimain_SOURCES = $(SRC) rcpi1main.f90
rpqdimain_FCFLAGS = -I../$(SRCOMMDIR)
rpqdimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librcpi.la
librpqdi_la_SOURCES = $(SRC)
librpqdi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpqdi_la_includedir = $(includedir)/$(SRINCDIR)
librpqdi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rcpi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
