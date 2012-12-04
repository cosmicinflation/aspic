SRC = cnbisr.f90 cnbireheat.f90
MOD = cnbisr.$(FC_MODEXT) cnbireheat.$(FC_MODEXT)

check_PROGRAMS = cnbimain
cnbimain_SOURCES = $(SRC) cnbimain.f90
cnbimain_FCFLAGS = -I../$(SRCOMMDIR)
cnbimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libcnbi.la
libcnbi_la_SOURCES = $(SRC)
libcnbi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libcnbi_la_includedir = $(includedir)/$(SRINCDIR)
libcnbi_la_include_HEADERS = $(MOD)

man_MANS = aspic_cnbi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



