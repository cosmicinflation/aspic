SRC = paisr.f90 paireheat.f90
MOD = paisr.$(FC_MODEXT) paireheat.$(FC_MODEXT)

check_PROGRAMS = paimain
paimain_SOURCES = $(SRC) paimain.f90
paimain_FCFLAGS = -I../$(SRCOMMDIR)
paimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libpai.la
libpai_la_SOURCES = $(SRC)
libpai_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libpai_la_includedir = $(includedir)/$(SRINCDIR)
libpai_la_include_HEADERS = $(MOD)

man_MANS = aspic_pai.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
