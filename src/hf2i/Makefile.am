SRC = hf2isr.f90 hf2ireheat.f90
MOD = hf2isr.$(FC_MODEXT) hf2ireheat.$(FC_MODEXT)

check_PROGRAMS = hf2imain
hf2imain_SOURCES = $(SRC) hf2imain.f90
hf2imain_FCFLAGS = -I../$(SRCOMMDIR)
hf2imain_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libhf2i.la
libhf2i_la_SOURCES = $(SRC)
libhf2i_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libhf2i_la_includedir = $(includedir)/$(SRINCDIR)
libhf2i_la_include_HEADERS = $(MOD)

man_MANS = aspic_hf2i.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



