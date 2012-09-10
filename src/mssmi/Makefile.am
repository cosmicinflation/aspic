SRC = mssmisr.f90 mssmireheat.f90
MOD = mssmisr.$(FC_MODEXT) mssmireheat.$(FC_MODEXT)

check_PROGRAMS = mssmimain
mssmimain_SOURCES = $(SRC) mssmimain.f90
mssmimain_FCFLAGS = -I../$(SRCOMMDIR)
mssmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libmssmi.la
libmssmi_la_SOURCES = $(SRC)
libmssmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libmssmi_la_includedir = $(includedir)/$(SRINCDIR)
libmssmi_la_include_HEADERS = $(MOD)

man_MANS = aspic_mssmi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



