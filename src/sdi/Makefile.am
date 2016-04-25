SRC = sdisr.f90 sdireheat.f90
MOD = sdisr.$(FC_MODEXT) sdireheat.$(FC_MODEXT)

check_PROGRAMS = sdimain
sdimain_SOURCES = $(SRC) sdimain.f90
sdimain_FCFLAGS = -I../$(SRCOMMDIR)
sdimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsdi.la
libsdi_la_SOURCES = $(SRC)
libsdi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsdi_la_includedir = $(includedir)/$(SRINCDIR)
libsdi_la_include_HEADERS = $(MOD)

man_MANS = aspic_sdi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
