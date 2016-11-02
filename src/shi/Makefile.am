SRC = shisr.f90 shireheat.f90
MOD = shisr.$(FC_MODEXT) shireheat.$(FC_MODEXT)

check_PROGRAMS = shimain
shimain_SOURCES = $(SRC) shimain.f90
shimain_FCFLAGS = -I../$(SRCOMMDIR)
shimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libshi.la
libshi_la_SOURCES = $(SRC)
libshi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libshi_la_includedir = $(includedir)/$(SRINCDIR)
libshi_la_include_HEADERS = $(MOD)

man_MANS = aspic_shi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
