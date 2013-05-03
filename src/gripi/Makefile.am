SRC = gripicommon.f90 ripisr.f90 ripireheat.f90 gripisr.F90 gripireheat.f90
MOD = gripicommon.$(FC_MODEXT) ripisr.$(FC_MODEXT) ripireheat.$(FC_MODEXT) \
	 gripisr.$(FC_MODEXT) gripireheat.$(FC_MODEXT)

check_PROGRAMS = ripimain gripimain
ripimain_SOURCES = $(SRC) ripimain.f90
ripimain_FCFLAGS = -I../$(SRCOMMDIR)
ripimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
gripimain_SOURCES = $(SRC) gripimain.f90
gripimain_FCFLAGS = -I../$(SRCOMMDIR)
gripimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libgripi.la
libgripi_la_SOURCES = $(SRC)
libgripi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libgripi_la_includedir = $(includedir)/$(SRINCDIR)
libgripi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ripi.3 aspic_gripi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
