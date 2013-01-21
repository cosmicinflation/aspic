SRC = gmssmicommon.f90 mssmisr.f90 mssmireheat.f90 gmssmisr.f90 gmssmireheat.f90
MOD = gmssmicommon.$(FC_MODEXT) mssmisr.$(FC_MODEXT) mssmireheat.$(FC_MODEXT) \
	gmssmisr.$(FC_MODEXT) gmssmireheat.$(FC_MODEXT)

check_PROGRAMS = mssmimain gmssmimain
mssmimain_SOURCES = $(SRC) mssmimain.f90
mssmimain_FCFLAGS = -I../$(SRCOMMDIR)
mssmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
gmssmimain_SOURCES = $(SRC) gmssmimain.f90
gmssmimain_FCFLAGS = -I../$(SRCOMMDIR)
gmssmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libgmssmi.la
libgmssmi_la_SOURCES = $(SRC)
libgmssmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libgmssmi_la_includedir = $(includedir)/$(SRINCDIR)
libgmssmi_la_include_HEADERS = $(MOD)

man_MANS = aspic_mssmi.3 aspic_gmssmi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
