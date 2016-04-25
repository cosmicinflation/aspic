SRC = sbkisr.f90 sbkireheat.f90
MOD = sbkisr.$(FC_MODEXT) sbkireheat.$(FC_MODEXT)

check_PROGRAMS = sbkimain
sbkimain_SOURCES = $(SRC) sbkimain.f90
sbkimain_FCFLAGS = -I../$(SRCOMMDIR)
sbkimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsbki.la
libsbki_la_SOURCES = $(SRC)
libsbki_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsbki_la_includedir = $(includedir)/$(SRINCDIR)
libsbki_la_include_HEADERS = $(MOD)

man_MANS = aspic_sbki.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
