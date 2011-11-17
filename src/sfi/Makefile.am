SRC = sfisr.f90 sfireheat.f90
MOD = sfisr.$(FC_MODEXT) sfireheat.$(FC_MODEXT)

check_PROGRAMS = sfimain
sfimain_SOURCES = $(SRC) sfimain.f90
sfimain_FCFLAGS = -I../$(SRCOMMDIR)
sfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsfi.la
libsfi_la_SOURCES = $(SRC)
libsfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsfi_la_includedir = $(includedir)/$(SRINCDIR)
libsfi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)