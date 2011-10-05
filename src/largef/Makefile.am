check_PROGRAMS = lfrehmain

lfrehmain_SOURCES = lfsrevol.f90 lfreheat.f90 lfrehmain.f90
lfrehmain_FCFLAGS = -I../$(SRCOMMDIR)
lfrehmain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a


clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)