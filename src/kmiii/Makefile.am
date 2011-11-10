SRC = kmiiisr.f90 kmiiireheat.f90
MOD = kmiiisr.$(FC_MODEXT) kmiiireheat.$(FC_MODEXT)

check_PROGRAMS = kmiiimain
kmiiimain_SOURCES = $(SRC) kmiiimain.f90
kmiiimain_FCFLAGS = -I../$(SRCOMMDIR)
kmiiimain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libkmiii.la
libkmiii_la_SOURCES = $(SRC)
libkmiii_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libkmiii_la_includedir = $(includedir)/$(SRINCDIR)
libkmiii_la_include_HEADERS = $(MOD)

clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)




