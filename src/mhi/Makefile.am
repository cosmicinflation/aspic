SRC = mhisr.f90 mhireheat.f90
MOD = mhisr.$(FC_MODEXT) mhireheat.$(FC_MODEXT)

check_PROGRAMS = mhimain
mhimain_SOURCES = $(SRC) mhimain.f90
mhimain_FCFLAGS = -I../$(SRCOMMDIR)
mhimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libmhi.la
libmhi_la_SOURCES = $(SRC)
libmhi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libmhi_la_includedir = $(includedir)/$(SRINCDIR)
libmhi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



