SRC = rchisr.f90 rchireheat.f90
MOD = rchisr.$(FC_MODEXT) rchireheat.$(FC_MODEXT)

check_PROGRAMS = rchimain
rchimain_SOURCES = $(SRC) rchimain.f90
rchimain_FCFLAGS = -I../$(SRCOMMDIR)
rchimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librchi.la
librchi_la_SOURCES = $(SRC)
librchi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librchi_la_includedir = $(includedir)/$(SRINCDIR)
librchi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rchi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
