SRC = kmiisr.f90 kmiireheat.f90
MOD = kmiisr.$(FC_MODEXT) kmiireheat.$(FC_MODEXT)

check_PROGRAMS = kmiimain
kmiimain_SOURCES = $(SRC) kmiimain.f90
kmiimain_FCFLAGS = -I../$(SRCOMMDIR)
kmiimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libkmii.la
libkmii_la_SOURCES = $(SRC)
libkmii_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libkmii_la_includedir = $(includedir)/$(SRINCDIR)
libkmii_la_include_HEADERS = $(MOD)

man_MANS = aspic_kmii.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



.NOTPARALLEL:
