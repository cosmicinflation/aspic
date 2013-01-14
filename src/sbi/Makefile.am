SRC = sbisr.f90 sbireheat.f90
MOD = sbisr.$(FC_MODEXT) sbireheat.$(FC_MODEXT)

check_PROGRAMS = sbimain
sbimain_SOURCES = $(SRC) sbimain.f90
sbimain_FCFLAGS = -I../$(SRCOMMDIR)
sbimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libsbi.la
libsbi_la_SOURCES = $(SRC)
libsbi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libsbi_la_includedir = $(includedir)/$(SRINCDIR)
libsbi_la_include_HEADERS = $(MOD)

man_MANS = aspic_sbi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



