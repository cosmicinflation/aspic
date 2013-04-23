SRC = wrhisr.f90 wrhireheat.f90
MOD = wrhisr.$(FC_MODEXT) wrhireheat.$(FC_MODEXT)

check_PROGRAMS = wrhimain
wrhimain_SOURCES = $(SRC) wrhimain.f90
wrhimain_FCFLAGS = -I../$(SRCOMMDIR)
wrhimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libwrhi.la
libwrhi_la_SOURCES = $(SRC)
libwrhi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libwrhi_la_includedir = $(includedir)/$(SRINCDIR)
libwrhi_la_include_HEADERS = $(MOD)

man_MANS = aspic_wrhi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
