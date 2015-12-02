SRC = hnisr.f90 hnireheat.f90
MOD = hnisr.$(FC_MODEXT) hnireheat.$(FC_MODEXT)

check_PROGRAMS = hnimain
hnimain_SOURCES = $(SRC) hnimain.f90
hnimain_FCFLAGS = -I../$(SRCOMMDIR)
hnimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhni.la
libhni_la_SOURCES = $(SRC)
libhni_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhni_la_includedir = $(includedir)/$(SRINCDIR)
libhni_la_include_HEADERS = $(MOD)

man_MANS = aspic_hni.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
