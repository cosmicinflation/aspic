SRC = ssicommon.f90 ssi1sr.f90 ssi1reheat.f90 \
                    ssi2sr.f90 ssi2reheat.f90 \
                    ssi3sr.f90 ssi3reheat.f90 \
                    ssi4sr.f90 ssi4reheat.f90 \
                    ssi5sr.f90 ssi5reheat.f90 \
                    ssi6sr.f90 ssi6reheat.f90 
MOD = ssicommon.$(FC_MODEXT)  ssi1sr.$(FC_MODEXT) ssi1reheat.$(FC_MODEXT) \
                              ssi2sr.$(FC_MODEXT) ssi2reheat.$(FC_MODEXT) \
                              ssi3sr.$(FC_MODEXT) ssi3reheat.$(FC_MODEXT) \
                              ssi4sr.$(FC_MODEXT) ssi4reheat.$(FC_MODEXT) \
                              ssi5sr.$(FC_MODEXT) ssi5reheat.$(FC_MODEXT) \
                              ssi6sr.$(FC_MODEXT) ssi6reheat.$(FC_MODEXT)

check_PROGRAMS = ssi1main ssi2main ssi3main ssi4main ssi5main ssi6main
ssi1main_SOURCES = $(SRC) ssi1main.f90
ssi1main_FCFLAGS = -I../$(SRCOMMDIR)
ssi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssi2main_SOURCES = $(SRC) ssi2main.f90
ssi2main_FCFLAGS = -I../$(SRCOMMDIR)
ssi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssi3main_SOURCES = $(SRC) ssi3main.f90
ssi3main_FCFLAGS = -I../$(SRCOMMDIR)
ssi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssi4main_SOURCES = $(SRC) ssi4main.f90
ssi4main_FCFLAGS = -I../$(SRCOMMDIR)
ssi4main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssi5main_SOURCES = $(SRC) ssi5main.f90
ssi5main_FCFLAGS = -I../$(SRCOMMDIR)
ssi5main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
ssi6main_SOURCES = $(SRC) ssi6main.f90
ssi6main_FCFLAGS = -I../$(SRCOMMDIR)
ssi6main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = libssi.la
libssi_la_SOURCES = $(SRC)
libssi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
libssi_la_includedir = $(includedir)/$(SRINCDIR)
libssi_la_include_HEADERS = $(MOD)

man_MANS = aspic_ssi1.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



