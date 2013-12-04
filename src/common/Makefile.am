COMMONSRC = infprec.F90 cosmopar.f90 infinout.f90 inftools.f90 \
            specialinf.F90 srflow.f90 srreheat.f90

COMMONMOD = cosmopar.$(FC_MODEXT) \
            infinout.$(FC_MODEXT) infprec.$(FC_MODEXT) inftools.$(FC_MODEXT) \
            hyp_2f1_module.$(FC_MODEXT) specialinf.$(FC_MODEXT) \
            specialprec.$(FC_MODEXT) srflow.$(FC_MODEXT) srreheat.$(FC_MODEXT)

check_LIBRARIES = libsrcommon.a
libsrcommon_a_SOURCES = $(COMMONSRC)

noinst_LTLIBRARIES = libsrcommon.la
libsrcommon_la_SOURCES = $(COMMONSRC)
#avoids race condition on names
libsrcommon_la_FCFLAGS = $(AM_FCFLAGS)

libsrcommon_la_includedir = $(includedir)/$(SRINCDIR)
libsrcommon_la_include_HEADERS = $(COMMONMOD)

#make clean does not remove modules by default
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)

.NOTPARALLEL:

