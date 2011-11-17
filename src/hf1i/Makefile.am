SRC = hf1isr.f90 hf1ireheat.f90
MOD = hf1isr.$(FC_MODEXT) hf1ireheat.$(FC_MODEXT)

check_PROGRAMS = hf1imain
hf1imain_SOURCES = $(SRC) hf1imain.f90
hf1imain_FCFLAGS = -I../$(SRCOMMDIR)
hf1imain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhf1i.la
libhf1i_la_SOURCES = $(SRC)
libhf1i_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhf1i_la_includedir = $(includedir)/$(SRINCDIR)
libhf1i_la_include_HEADERS = $(MOD)

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



