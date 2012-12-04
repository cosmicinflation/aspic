SRC = ripisr.f90 ripireheat.f90
MOD = ripisr.$(FC_MODEXT) ripireheat.$(FC_MODEXT)

check_PROGRAMS = ripimain
ripimain_SOURCES = $(SRC) ripimain.f90
ripimain_FCFLAGS = -I../$(SRCOMMDIR)
ripimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libripi.la
libripi_la_SOURCES = $(SRC)
libripi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libripi_la_includedir = $(includedir)/$(SRINCDIR)
libripi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ripi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



