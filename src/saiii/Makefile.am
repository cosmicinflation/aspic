SRC = saiiicommon.f90 saiiicomreh.f90 saiii1sr.f90 saiii1reheat.f90 \
	saiii2sr.f90 saiii2reheat.f90 saiii3sr.f90 saiii3reheat.f90
MOD = saiiicommon.$(FC_MODEXT) saiiicomreh.$(FC_MODEXT) \
	saiii1sr.$(FC_MODEXT) saiii1reheat.$(FC_MODEXT) \
	saiii2sr.$(FC_MODEXT) saiii2reheat.$(FC_MODEXT) \
	saiii3sr.$(FC_MODEXT) saiii3reheat.$(FC_MODEXT)

check_PROGRAMS = saiii1main saiii2main saiii3main
saiii1main_SOURCES = $(SRC) saiii1main.f90
saiii1main_FCFLAGS = -I../$(SRCOMMDIR)
saiii1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
saiii2main_SOURCES = $(SRC) saiii2main.f90
saiii2main_FCFLAGS = -I../$(SRCOMMDIR)
saiii2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
saiii3main_SOURCES = $(SRC) saiii3main.f90
saiii3main_FCFLAGS = -I../$(SRCOMMDIR)
saiii3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsaiii.la
libsaiii_la_SOURCES = $(SRC)
libsaiii_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsaiii_la_includedir = $(includedir)/$(SRINCDIR)
libsaiii_la_include_HEADERS = $(MOD)

man_MANS = aspic_saiii1.3 aspic_saiii2.3 aspic_saiii3.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
