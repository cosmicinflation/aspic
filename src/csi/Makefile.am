SRC = csisr.f90 csireheat.f90
MOD = csisr.$(FC_MODEXT) csireheat.$(FC_MODEXT)

check_PROGRAMS = csimain
csimain_SOURCES = $(SRC) csimain.f90
csimain_FCFLAGS = -I../$(SRCOMMDIR)
csimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libcsi.la
libcsi_la_SOURCES = $(SRC)
libcsi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libcsi_la_includedir = $(includedir)/$(SRINCDIR)
libcsi_la_include_HEADERS = $(MOD)

man_MANS = aspic_csi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
