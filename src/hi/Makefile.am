SRC = hisr.f90 hireheat.f90
MOD = hisr.$(FC_MODEXT) hireheat.$(FC_MODEXT)

check_PROGRAMS = himain
himain_SOURCES = $(SRC) himain.f90
himain_FCFLAGS = -I../$(SRCOMMDIR)
himain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhi.la
libhi_la_SOURCES = $(SRC)
libhi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhi_la_includedir = $(includedir)/$(SRINCDIR)
libhi_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



