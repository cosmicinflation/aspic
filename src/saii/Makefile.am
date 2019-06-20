SRC = saiicommon.f90 saiicomreh.f90 saii1sr.f90 saii1reheat.f90 \
	saii2sr.f90 saii2reheat.f90

MOD = saiicommon.$(FC_MODEXT) saiicomreh.$(FC_MODEXT) \
	saii1sr.$(FC_MODEXT) saii1reheat.$(FC_MODEXT) \
	saii2sr.$(FC_MODEXT) saii2reheat.$(FC_MODEXT) 

check_PROGRAMS = saii1main
saii1main_SOURCES = $(SRC) saii1main.f90
saii1main_FCFLAGS = -I../$(SRCOMMDIR)
saii1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

check_PROGRAMS = saii2main
saii2main_SOURCES = $(SRC) saii2main.f90
saii2main_FCFLAGS = -I../$(SRCOMMDIR)
saii2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsaii.la
libsaii_la_SOURCES = $(SRC)
libsaii_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsaii_la_includedir = $(includedir)/$(SRINCDIR)
libsaii_la_include_HEADERS = $(MOD)

man_MANS = aspic_saii1.3 aspic_saii2.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
