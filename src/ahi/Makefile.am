SRC = ahisr.f90 ahireheat.f90
MOD = ahisr.$(FC_MODEXT) ahireheat.$(FC_MODEXT)

check_PROGRAMS = ahimain
ahimain_SOURCES = $(SRC) ahimain.f90
ahimain_FCFLAGS = -I../$(SRCOMMDIR)
ahimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libahi.la
libahi_la_SOURCES = $(SRC)
libahi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libahi_la_includedir = $(includedir)/$(SRINCDIR)
libahi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ahi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
