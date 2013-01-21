SRC = rcqisr.f90 rcqireheat.f90
MOD = rcqisr.$(FC_MODEXT) rcqireheat.$(FC_MODEXT)

check_PROGRAMS = rcqimain
rcqimain_SOURCES = $(SRC) rcqimain.f90
rcqimain_FCFLAGS = -I../$(SRCOMMDIR)
rcqimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librcqi.la
librcqi_la_SOURCES = $(SRC)
librcqi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librcqi_la_includedir = $(includedir)/$(SRINCDIR)
librcqi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rcqi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
