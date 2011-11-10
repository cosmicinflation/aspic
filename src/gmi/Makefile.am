SRC = gmisr.f90 gmireheat.f90
MOD = gmisr.$(FC_MODEXT) gmireheat.$(FC_MODEXT)

check_PROGRAMS = gmimain
gmimain_SOURCES = $(SRC) gmimain.f90
gmimain_FCFLAGS = -I../$(SRCOMMDIR)
gmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libgmi.la
libgmi_la_SOURCES = $(SRC)
libgmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libgmi_la_includedir = $(includedir)/$(SRINCDIR)
libgmi_la_include_HEADERS = $(MOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)




