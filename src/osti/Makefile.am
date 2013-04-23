SRC = ostisr.f90 ostireheat.f90
MOD = ostisr.$(FC_MODEXT) ostireheat.$(FC_MODEXT)

check_PROGRAMS = ostimain
ostimain_SOURCES = $(SRC) ostimain.f90
ostimain_FCFLAGS = -I../$(SRCOMMDIR)
ostimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libosti.la
libosti_la_SOURCES = $(SRC)
libosti_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libosti_la_includedir = $(includedir)/$(SRINCDIR)
libosti_la_include_HEADERS = $(MOD)

man_MANS = aspic_osti.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
