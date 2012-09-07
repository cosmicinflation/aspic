SRC = rmisr.f90 rmireheat.f90
MOD = rmisr.$(FC_MODEXT) rmireheat.$(FC_MODEXT)

check_PROGRAMS = rmimain
rmimain_SOURCES = $(SRC) rmimain.f90
rmimain_FCFLAGS = -I../$(SRCOMMDIR)
rmimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librmi.la
librmi_la_SOURCES = $(SRC)
librmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librmi_la_includedir = $(includedir)/$(SRINCDIR)
librmi_la_include_HEADERS = $(MOD)

man_MANS = librmi.3

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
