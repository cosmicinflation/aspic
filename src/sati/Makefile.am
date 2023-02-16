SRC = satisr.f90 satireheat.f90
MOD = satisr.$(FC_MODEXT) satireheat.$(FC_MODEXT)

check_PROGRAMS = satimain
satimain_SOURCES = $(SRC) satimain.f90
satimain_FCFLAGS = -I../$(SRCOMMDIR)
satimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsati.la
libsati_la_SOURCES = $(SRC)
libsati_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsati_la_includedir = $(includedir)/$(SRINCDIR)
libsati_la_include_HEADERS = $(MOD)

man_MANS = aspic_sati.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
