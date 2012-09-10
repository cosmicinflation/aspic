SRC = nckisr.f90 nckireheat.f90
MOD = nckisr.$(FC_MODEXT) nckireheat.$(FC_MODEXT)

check_PROGRAMS = nckimain
nckimain_SOURCES = $(SRC) nckimain.f90
nckimain_FCFLAGS = -I../$(SRCOMMDIR)
nckimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libncki.la
libncki_la_SOURCES = $(SRC)
libncki_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libncki_la_includedir = $(includedir)/$(SRINCDIR)
libncki_la_include_HEADERS = $(MOD)

man_MANS = aspic_ncki.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



