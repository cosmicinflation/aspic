SRC = lisr.f90 lireheat.f90
MOD = lisr.$(FC_MODEXT) lireheat.$(FC_MODEXT)

check_PROGRAMS = limain
limain_SOURCES = $(SRC) limain.f90
limain_FCFLAGS = -I../$(SRCOMMDIR)
limain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libli.la
libli_la_SOURCES = $(SRC)
libli_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libli_la_includedir = $(includedir)/$(SRINCDIR)
libli_la_include_HEADERS = $(MOD)

man_MANS = aspic_li.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



