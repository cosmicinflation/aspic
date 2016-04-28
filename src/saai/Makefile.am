SRC = saaisr.f90 saaireheat.f90
MOD = saaisr.$(FC_MODEXT) saaireheat.$(FC_MODEXT)

check_PROGRAMS = saaimain
saaimain_SOURCES = $(SRC) saaimain.f90
saaimain_FCFLAGS = -I../$(SRCOMMDIR)
saaimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsaai.la
libsaai_la_SOURCES = $(SRC)
libsaai_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsaai_la_includedir = $(includedir)/$(SRINCDIR)
libsaai_la_include_HEADERS = $(MOD)

man_MANS = aspic_saai.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
