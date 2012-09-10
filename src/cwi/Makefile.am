SRC = cwisr.f90 cwireheat.f90
MOD = cwisr.$(FC_MODEXT) cwireheat.$(FC_MODEXT)

check_PROGRAMS = cwimain
cwimain_SOURCES = $(SRC) cwimain.f90
cwimain_FCFLAGS = -I../$(SRCOMMDIR)
cwimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libcwi.la
libcwi_la_SOURCES = $(SRC)
libcwi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libcwi_la_includedir = $(includedir)/$(SRINCDIR)
libcwi_la_include_HEADERS = $(MOD)

man_MANS = aspic_cwi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



