SRC = deisr.f90 deireheat.f90
MOD = deisr.$(FC_MODEXT) deireheat.$(FC_MODEXT)

check_PROGRAMS = deimain
deimain_SOURCES = $(SRC) deimain.f90
deimain_FCFLAGS = -I../$(SRCOMMDIR)
deimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libdei.la
libdei_la_SOURCES = $(SRC)
libdei_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libdei_la_includedir = $(includedir)/$(SRINCDIR)
libdei_la_include_HEADERS = $(MOD)

man_MANS = aspic_dei.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
