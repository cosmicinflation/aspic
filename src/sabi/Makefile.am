SRC = sabisr.f90 sabireheat.f90
MOD = sabisr.$(FC_MODEXT) sabireheat.$(FC_MODEXT)

check_PROGRAMS = sabimain
sabimain_SOURCES = $(SRC) sabimain.f90
sabimain_FCFLAGS = -I../$(SRCOMMDIR)
sabimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsabi.la
libsabi_la_SOURCES = $(SRC)
libsabi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsabi_la_includedir = $(includedir)/$(SRINCDIR)
libsabi_la_include_HEADERS = $(MOD)

man_MANS = aspic_sabi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)
.NOTPARALLEL:
