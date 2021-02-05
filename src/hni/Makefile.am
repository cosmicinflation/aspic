SRC = hnicommon.f90 hnicomreh.f90 hni1sr.f90 hni1reheat.f90 \
	hni2sr.f90 hni2reheat.f90

MOD = hnicommon.$(FC_MODEXT) hnicomreh.$(FC_MODEXT) \
	hni1sr.$(FC_MODEXT) hni1reheat.$(FC_MODEXT) \
	hni2sr.$(FC_MODEXT) hni2reheat.$(FC_MODEXT)

check_PROGRAMS = hni1main hni2main
hni1main_SOURCES = $(SRC) hni1main.f90
hni1main_FCFLAGS = -I../$(SRCOMMDIR)
hni1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
hni2main_SOURCES = $(SRC) hni2main.f90
hni2main_FCFLAGS = -I../$(SRCOMMDIR)
hni2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhni.la
libhni_la_SOURCES = $(SRC)
libhni_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhni_la_includedir = $(includedir)/$(SRINCDIR)
libhni_la_include_HEADERS = $(MOD)

man_MANS = aspic_hni1.3 aspic_hni2.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
