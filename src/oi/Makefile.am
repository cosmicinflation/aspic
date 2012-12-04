SRC = oisr.f90 oireheat.f90
MOD = oisr.$(FC_MODEXT) oireheat.$(FC_MODEXT)

check_PROGRAMS = oimain
oimain_SOURCES = $(SRC) oimain.f90
oimain_FCFLAGS = -I../$(SRCOMMDIR)
oimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liboi.la
liboi_la_SOURCES = $(SRC)
liboi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liboi_la_includedir = $(includedir)/$(SRINCDIR)
liboi_la_include_HEADERS = $(MOD)

man_MANS = aspic_oi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
