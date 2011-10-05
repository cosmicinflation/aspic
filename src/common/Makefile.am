COMMONSRC = infprec.f90 cosmopar.f90 infinout.f90 inftools.f90 srreheat.f90

check_LIBRARIES = libsrcommon.a
libsrcommon_a_SOURCES = $(COMMONSRC)

noinst_LTLIBRARIES = libsrcommon.la
libsrcommon_la_SOURCES = $(COMMONSRC)
#avoids race condition on names
libsrcommon_la_FCFLAGS = $(AM_FCFLAGS)

#make clean does not remove modules by default
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)