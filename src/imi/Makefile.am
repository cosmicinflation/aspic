SRC = imisr.f90 imireheat.f90
MOD = imisr.$(FC_MODEXT) imireheat.$(FC_MODEXT)

check_PROGRAMS = imimain
imimain_SOURCES = $(SRC) imimain.f90
imimain_FCFLAGS = -I../$(SRCOMMDIR)
imimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libimi.la
libimi_la_SOURCES = $(SRC)
libimi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libimi_la_includedir = $(includedir)/$(SRINCDIR)
libimi_la_include_HEADERS = $(MOD)

man_MANS = aspic_imi.3


clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)

.NOTPARALLEL:
