SRC = oripisr.f90 oripireheat.f90
MOD = oripisr.$(FC_MODEXT) oripireheat.$(FC_MODEXT)

check_PROGRAMS = oripimain
oripimain_SOURCES = $(SRC) oripimain.f90
oripimain_FCFLAGS = -I../$(SRCOMMDIR)
oripimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = liboripi.la
liboripi_la_SOURCES = $(SRC)
liboripi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
liboripi_la_includedir = $(includedir)/$(SRINCDIR)
liboripi_la_include_HEADERS = $(MOD)

man_MANS = aspic_oripi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
