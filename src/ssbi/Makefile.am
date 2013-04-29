SRC = ssbicommon.f90 ssbicomreh.f90 \
	ssbi1sr.f90 ssbi1reheat.f90 \
	ssbi2sr.f90 ssbi2reheat.f90 \
	ssbi3sr.f90 ssbi3reheat.f90 \
	ssbi4sr.f90 ssbi4reheat.f90 \
	ssbi5sr.f90 ssbi5reheat.f90 \
	ssbi6sr.f90 ssbi6reheat.f90 
MOD = ssbicommon.$(FC_MODEXT) ssbicomreh.$(FC_MODEXT) \
	ssbi1sr.$(FC_MODEXT) ssbi1reheat.$(FC_MODEXT) \
	ssbi2sr.$(FC_MODEXT) ssbi2reheat.$(FC_MODEXT) \
	ssbi3sr.$(FC_MODEXT) ssbi3reheat.$(FC_MODEXT) \
	ssbi4sr.$(FC_MODEXT) ssbi4reheat.$(FC_MODEXT) \
	ssbi5sr.$(FC_MODEXT) ssbi5reheat.$(FC_MODEXT) \
	ssbi6sr.$(FC_MODEXT) ssbi6reheat.$(FC_MODEXT)

check_PROGRAMS = ssbi1main ssbi2main ssbi3main ssbi4main ssbi5main ssbi6main
ssbi1main_SOURCES = $(SRC) ssbi1main.f90
ssbi1main_FCFLAGS = -I../$(SRCOMMDIR)
ssbi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssbi2main_SOURCES = $(SRC) ssbi2main.f90
ssbi2main_FCFLAGS = -I../$(SRCOMMDIR)
ssbi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssbi3main_SOURCES = $(SRC) ssbi3main.f90
ssbi3main_FCFLAGS = -I../$(SRCOMMDIR)
ssbi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssbi4main_SOURCES = $(SRC) ssbi4main.f90
ssbi4main_FCFLAGS = -I../$(SRCOMMDIR)
ssbi4main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssbi5main_SOURCES = $(SRC) ssbi5main.f90
ssbi5main_FCFLAGS = -I../$(SRCOMMDIR)
ssbi5main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssbi6main_SOURCES = $(SRC) ssbi6main.f90
ssbi6main_FCFLAGS = -I../$(SRCOMMDIR)
ssbi6main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libssbi.la
libssbi_la_SOURCES = $(SRC)
libssbi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libssbi_la_includedir = $(includedir)/$(SRINCDIR)
libssbi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ssbi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
