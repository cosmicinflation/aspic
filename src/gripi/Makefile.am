SRC = gripisr.f90 gripireheat.f90
MOD = gripisr.$(FC_MODEXT) gripireheat.$(FC_MODEXT)

check_PROGRAMS = gripimain
gripimain_SOURCES = $(SRC) gripimain.f90
gripimain_FCFLAGS = -I../$(SRCOMMDIR)
gripimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libgripi.la
libgripi_la_SOURCES = $(SRC)
libgripi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libgripi_la_includedir = $(includedir)/$(SRINCDIR)
libgripi_la_include_HEADERS = $(MOD)

man_MANS = aspic_gripi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
