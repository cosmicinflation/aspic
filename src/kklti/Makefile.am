SRC = kkltisr.f90 kkltireheat.f90
MOD = kkltisr.$(FC_MODEXT) kkltireheat.$(FC_MODEXT)

check_PROGRAMS = kkltimain
kkltimain_SOURCES = $(SRC) kkltimain.f90
kkltimain_FCFLAGS = -I../$(SRCOMMDIR)
kkltimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libkklti.la
libkklti_la_SOURCES = $(SRC)
libkklti_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libkklti_la_includedir = $(includedir)/$(SRINCDIR)
libkklti_la_include_HEADERS = $(MOD)

man_MANS = aspic_kklti.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
