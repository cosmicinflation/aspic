SRC = nisr.f90 nireheat.f90
MOD = nisr.$(FC_MODEXT) nireheat.$(FC_MODEXT)

check_PROGRAMS = nimain
nimain_SOURCES = $(SRC) nimain.f90
nimain_FCFLAGS = -I../$(SRCOMMDIR)
nimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libni.la
libni_la_SOURCES = $(SRC)
libni_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libni_la_includedir = $(includedir)/$(SRINCDIR)
libni_la_include_HEADERS = $(MOD)

man_MANS = aspic_ni.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
