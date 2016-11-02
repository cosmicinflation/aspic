SRC = rpqtisr.f90 rpqtireheat.f90
MOD = rpqtisr.$(FC_MODEXT) rpqtireheat.$(FC_MODEXT)

check_PROGRAMS = rpqtimain
rpqtimain_SOURCES = $(SRC) rpqtimain.f90
rpqtimain_FCFLAGS = -I../$(SRCOMMDIR)
rpqtimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librpqti.la
librpqti_la_SOURCES = $(SRC)
librpqti_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librpqti_la_includedir = $(includedir)/$(SRINCDIR)
librpqti_la_include_HEADERS = $(MOD)

man_MANS = aspic_rpqti.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
