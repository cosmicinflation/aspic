SRC = cnaisr.f90 cnaireheat.f90
MOD = cnaisr.$(FC_MODEXT) cnaireheat.$(FC_MODEXT)

check_PROGRAMS = cnaimain
cnaimain_SOURCES = $(SRC) cnaimain.f90
cnaimain_FCFLAGS = -I../$(SRCOMMDIR)
cnaimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libcnai.la
libcnai_la_SOURCES = $(SRC)
libcnai_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libcnai_la_includedir = $(includedir)/$(SRINCDIR)
libcnai_la_include_HEADERS = $(MOD)

man_MANS = aspic_cnai.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
