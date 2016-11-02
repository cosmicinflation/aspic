SRC = rpqdisr.f90 rpqdireheat.f90
MOD = rpqdisr.$(FC_MODEXT) rpqdireheat.$(FC_MODEXT)

check_PROGRAMS = rpqdimain
rpqdimain_SOURCES = $(SRC) rpqdimain.f90
rpqdimain_FCFLAGS = -I../$(SRCOMMDIR)
rpqdimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librpqdi.la
librpqdi_la_SOURCES = $(SRC)
librpqdi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpqdi_la_includedir = $(includedir)/$(SRINCDIR)
librpqdi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rpqdi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
