SRC = aisr.f90 aireheat.f90
MOD = aisr.$(FC_MODEXT) aireheat.$(FC_MODEXT)

check_PROGRAMS = aimain
aimain_SOURCES = $(SRC) aimain.f90
aimain_FCFLAGS = -I../$(SRCOMMDIR)
aimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libai.la
libai_la_SOURCES = $(SRC)
libai_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libai_la_includedir = $(includedir)/$(SRINCDIR)
libai_la_include_HEADERS = $(MOD)

man_MANS = aspic_ai.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



