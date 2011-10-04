check_PROGRAMS = sfrehmain

sfrehmain_SOURCES = sfsrevol.f90 sfreheat.f90 sfrehmain.f90
sfrehmain_FCFLAGS = -I../$(SRCOMMDIR)
sfrehmain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a


clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)