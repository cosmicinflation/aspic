SRC = bsusybisr.f90 bsusybireheat.f90
MOD = bsusybisr.$(FC_MODEXT) bsusybireheat.$(FC_MODEXT)

check_PROGRAMS = bsusybimain
bsusybimain_SOURCES = $(SRC) bsusybimain.f90
bsusybimain_FCFLAGS = -I../$(SRCOMMDIR)
bsusybimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libbsusybi.la
libbsusybi_la_SOURCES = $(SRC)
libbsusybi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libbsusybi_la_includedir = $(includedir)/$(SRINCDIR)
libbsusybi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
