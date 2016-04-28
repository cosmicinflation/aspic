SRC = sacisr.f90 sacireheat.f90
MOD = sacisr.$(FC_MODEXT) sacireheat.$(FC_MODEXT)

check_PROGRAMS = sacimain
sacimain_SOURCES = $(SRC) sacimain.f90
sacimain_FCFLAGS = -I../$(SRCOMMDIR)
sacimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsaci.la
libsaci_la_SOURCES = $(SRC)
libsaci_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsaci_la_includedir = $(includedir)/$(SRINCDIR)
libsaci_la_include_HEADERS = $(MOD)

man_MANS = aspic_saci.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
