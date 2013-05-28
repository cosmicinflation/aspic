SRC = wrisr.f90 wrireheat.f90
MOD = wrisr.$(FC_MODEXT) wrireheat.$(FC_MODEXT)

check_PROGRAMS = wrimain
wrimain_SOURCES = $(SRC) wrimain.f90
wrimain_FCFLAGS = -I../$(SRCOMMDIR)
wrimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libwri.la
libwri_la_SOURCES = $(SRC)
libwri_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libwri_la_includedir = $(includedir)/$(SRINCDIR)
libwri_la_include_HEADERS = $(MOD)

man_MANS = aspic_wri.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
