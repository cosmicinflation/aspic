SRC = beisr.f90 beireheat.f90
MOD = beisr.$(FC_MODEXT) beireheat.$(FC_MODEXT)

check_PROGRAMS = beimain
beimain_SOURCES = $(SRC) beimain.f90
beimain_FCFLAGS = -I../$(SRCOMMDIR)
beimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libbei.la
libbei_la_SOURCES = $(SRC)
libbei_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libbei_la_includedir = $(includedir)/$(SRINCDIR)
libbei_la_include_HEADERS = $(MOD)

man_MANS = aspic_bei.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
