SRC = rcmisr.f90 rcmireheat.f90
MOD = rcmisr.$(FC_MODEXT) rcmireheat.$(FC_MODEXT)

check_PROGRAMS = rcmimain
rcmimain_SOURCES = $(SRC) rcmimain.f90
rcmimain_FCFLAGS = -I../$(SRCOMMDIR)
rcmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librcmi.la
librcmi_la_SOURCES = $(SRC)
librcmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librcmi_la_includedir = $(includedir)/$(SRINCDIR)
librcmi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



