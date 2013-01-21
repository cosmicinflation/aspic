SRC = dsisr.f90 dsireheat.f90
MOD = dsisr.$(FC_MODEXT) dsireheat.$(FC_MODEXT)

check_PROGRAMS = dsimain
dsimain_SOURCES = $(SRC) dsimain.f90
dsimain_FCFLAGS = -I../$(SRCOMMDIR)
dsimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libdsi.la
libdsi_la_SOURCES = $(SRC)
libdsi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libdsi_la_includedir = $(includedir)/$(SRINCDIR)
libdsi_la_include_HEADERS = $(MOD)

man_MANS = aspic_dsi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
