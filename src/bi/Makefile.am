SRC = bisr.f90 bireheat.f90
MOD =  bisr.$(FC_MODEXT) bireheat.$(FC_MODEXT)

check_PROGRAMS = bimain
bimain_SOURCES = $(SRC) bimain.f90
bimain_FCFLAGS = -I../$(SRCOMMDIR)
bimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libbi.la
libbi_la_SOURCES = $(SRC)
libbi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libbi_la_includedir = $(includedir)/$(SRINCDIR)
libbi_la_include_HEADERS = $(MOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
