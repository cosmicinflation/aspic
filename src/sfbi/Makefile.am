SRC = sfbisr.f90 sfbireheat.f90
MOD = sfbisr.$(FC_MODEXT) sfbireheat.$(FC_MODEXT)

check_PROGRAMS = sfbimain
sfbimain_SOURCES = $(SRC) sfbimain.f90
sfbimain_FCFLAGS = -I../$(SRCOMMDIR) -I../sfi
sfbimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a 
sfbimain_LDADD += ../sfi/libsfi_la-sfireheat.o
sfbimain_LDADD += ../sfi/libsfi_la-sfisr.o

noinst_LTLIBRARIES = libsfbi.la
libsfbi_la_SOURCES = $(SRC)
libsfbi_la_FCFLAGS = -I../$(SRCOMMDIR) -I../sfi $(AM_FCFLAGS)
libsfbi_la_includedir = $(includedir)/$(SRINCDIR)
libsfbi_la_include_HEADERS = $(MOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
