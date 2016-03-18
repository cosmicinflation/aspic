SRC = vfmieos.f90 vfmisr.f90 vfmireheat.f90

MOD = vfmieos.$(FC_MODEXT) vfmisr.$(FC_MODEXT) vfmireheat.$(FC_MODEXT)

check_PROGRAMS = vfmimain
vfmimain_SOURCES = $(SRC) vfmimain.f90
vfmimain_FCFLAGS = -I../$(SRCOMMDIR)
vfmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libvfmi.la
libvfmi_la_SOURCES = $(SRC)
libvfmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libvfmi_la_includedir = $(includedir)/$(SRINCDIR)
libvfmi_la_include_HEADERS = $(MOD)

man_MANS = aspic_vfmi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
