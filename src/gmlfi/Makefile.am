SRC = gmlfisr.f90 gmlfireheat.f90
MOD = gmlfisr.$(FC_MODEXT) gmlfireheat.$(FC_MODEXT)

check_PROGRAMS = mlfimain gmlfimain
gmlfimain_SOURCES = $(SRC) gmlfimain.f90
gmlfimain_FCFLAGS = -I../$(SRCOMMDIR)
gmlfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
mlfimain_SOURCES = $(SRC) mlfimain.f90
mlfimain_FCFLAGS = -I../$(SRCOMMDIR)
mlfimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libgmlfi.la
libgmlfi_la_SOURCES = $(SRC)
libgmlfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libgmlfi_la_includedir = $(includedir)/$(SRINCDIR)
libgmlfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_gmlfi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
