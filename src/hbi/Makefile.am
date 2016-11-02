SRC = hbisr.f90 hbireheat.f90
MOD = hbisr.$(FC_MODEXT) hbireheat.$(FC_MODEXT)

check_PROGRAMS = hbimain
hbimain_SOURCES = $(SRC) hbimain.f90
hbimain_FCFLAGS = -I../$(SRCOMMDIR)
hbimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhbi.la
libhbi_la_SOURCES = $(SRC)
libhbi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhbi_la_includedir = $(includedir)/$(SRINCDIR)
libhbi_la_include_HEADERS = $(MOD)

man_MANS = aspic_hbi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
