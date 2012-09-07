SRC = mlfisr.f90 mlfireheat.f90
MOD = mlfisr.$(FC_MODEXT) mlfireheat.$(FC_MODEXT)

check_PROGRAMS = mlfimain
mlfimain_SOURCES = $(SRC) mlfimain.f90
mlfimain_FCFLAGS = -I../$(SRCOMMDIR)
mlfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libmlfi.la
libmlfi_la_SOURCES = $(SRC)
libmlfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libmlfi_la_includedir = $(includedir)/$(SRINCDIR)
libmlfi_la_include_HEADERS = $(MOD)

man_MANS = libmlfi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)