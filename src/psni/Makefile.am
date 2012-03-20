SRC = psnisr.f90 psnireheat.f90
MOD = psnisr.$(FC_MODEXT) psnireheat.$(FC_MODEXT)

check_PROGRAMS = psnimain
psnimain_SOURCES = $(SRC) psnimain.f90
psnimain_FCFLAGS = -I../$(SRCOMMDIR)
psnimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libpsni.la
libpsni_la_SOURCES = $(SRC)
libpsni_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libpsni_la_includedir = $(includedir)/$(SRINCDIR)
libpsni_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



