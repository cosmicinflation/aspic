SRC = nficommon.f90 nficomreh.f90 \
	nfi1sr.f90 nfi1reheat.f90 \
	nfi2sr.f90 nfi2reheat.f90 \
	nfi3sr.f90 nfi3reheat.f90 \
	nfi4sr.f90 nfi4reheat.f90

MOD = nficommon.$(FC_MODEXT) nficomreh.$(FC_MODEXT) \
	nfi1sr.$(FC_MODEXT) nfi1reheat.$(FC_MODEXT) \
	nfi2sr.$(FC_MODEXT) nfi2reheat.$(FC_MODEXT) \
	nfi3sr.$(FC_MODEXT) nfi3reheat.$(FC_MODEXT) \
	nfi4sr.$(FC_MODEXT) nfi4reheat.$(FC_MODEXT)


check_PROGRAMS = nfi1main nfi2main nfi3main nfi4main
nfi1main_SOURCES = $(SRC) nfi1main.f90
nfi1main_FCFLAGS = -I../$(SRCOMMDIR)
nfi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
nfi2main_SOURCES = $(SRC) nfi2main.f90
nfi2main_FCFLAGS = -I../$(SRCOMMDIR)
nfi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
nfi3main_SOURCES = $(SRC) nfi3main.f90
nfi3main_FCFLAGS = -I../$(SRCOMMDIR)
nfi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
nfi4main_SOURCES = $(SRC) nfi4main.f90
nfi4main_FCFLAGS = -I../$(SRCOMMDIR)
nfi4main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libnfi.la
libnfi_la_SOURCES = $(SRC)
libnfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libnfi_la_includedir = $(includedir)/$(SRINCDIR)
libnfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_nfi1.3 aspic_nfi2.3 aspic_nfi3.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
