noinst_LIBRARIES = libsrcommon.a
libsrcommon_a_SOURCES = infprec.f90 cosmopar.f90 infinout.f90 inftools.f90 srreheat.f90

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)