SRC = cndisr.f90 cndireheat.f90
MOD = cndisr.$(FC_MODEXT) cndireheat.$(FC_MODEXT)

check_PROGRAMS = cndimain
cndimain_SOURCES = $(SRC) cndimain.f90
cndimain_FCFLAGS = -I../$(SRCOMMDIR)
cndimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libcndi.la
libcndi_la_SOURCES = $(SRC)
libcndi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libcndi_la_includedir = $(includedir)/$(SRINCDIR)
libcndi_la_include_HEADERS = $(MOD)

man_MANS = aspic_cndi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



