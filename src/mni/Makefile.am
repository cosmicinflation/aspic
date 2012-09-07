SRC = mnisr.f90 mnireheat.f90
MOD = mnisr.$(FC_MODEXT) mnireheat.$(FC_MODEXT)

check_PROGRAMS = mnimain
mnimain_SOURCES = $(SRC) mnimain.f90
mnimain_FCFLAGS = -I../$(SRCOMMDIR)
mnimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libmni.la
libmni_la_SOURCES = $(SRC)
libmni_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libmni_la_includedir = $(includedir)/$(SRINCDIR)
libmni_la_include_HEADERS = $(MOD)

man_MANS = libmni.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



