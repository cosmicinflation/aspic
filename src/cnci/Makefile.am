SRC = cncisr.f90 cncireheat.f90
MOD = cncisr.$(FC_MODEXT) cncireheat.$(FC_MODEXT)

check_PROGRAMS = cncimain
cncimain_SOURCES = $(SRC) cncimain.f90
cncimain_FCFLAGS = -I../$(SRCOMMDIR)
cncimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libcnci.la
libcnci_la_SOURCES = $(SRC)
libcnci_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libcnci_la_includedir = $(includedir)/$(SRINCDIR)
libcnci_la_include_HEADERS = $(MOD)

man_MANS = aspic_cnci.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
