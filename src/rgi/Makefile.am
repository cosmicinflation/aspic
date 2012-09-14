SRC = rgisr.f90 rgireheat.f90
MOD = rgisr.$(FC_MODEXT) rgireheat.$(FC_MODEXT)

check_PROGRAMS = rgimain
rgimain_SOURCES = $(SRC) rgimain.f90
rgimain_FCFLAGS = -I../$(SRCOMMDIR)
rgimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librgi.la
librgi_la_SOURCES = $(SRC)
librgi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librgi_la_includedir = $(includedir)/$(SRINCDIR)
librgi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)


