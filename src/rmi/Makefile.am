SRC = rmicommon.f90 rmi1sr.f90 rmi1reheat.f90 \
                    rmi2sr.f90 rmi2reheat.f90 \
                    rmi3sr.f90 rmi3reheat.f90 \
                    rmi4sr.f90 rmi4reheat.f90 
MOD = rmicommon.$(FC_MODEXT)  rmi1sr.$(FC_MODEXT) rmi1reheat.$(FC_MODEXT) \
                              rmi2sr.$(FC_MODEXT) rmi2reheat.$(FC_MODEXT) \
                              rmi3sr.$(FC_MODEXT) rmi3reheat.$(FC_MODEXT) \
                              rmi4sr.$(FC_MODEXT) rmi4reheat.$(FC_MODEXT) 

check_PROGRAMS = rmi1main rmi2main rmi3main rmi4main
rmi1main_SOURCES = $(SRC) rmi1main.f90
rmi1main_FCFLAGS = -I../$(SRCOMMDIR)
rmi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
rmi2main_SOURCES = $(SRC) rmi2main.f90
rmi2main_FCFLAGS = -I../$(SRCOMMDIR)
rmi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
rmi3main_SOURCES = $(SRC) rmi3main.f90
rmi3main_FCFLAGS = -I../$(SRCOMMDIR)
rmi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
rmi4main_SOURCES = $(SRC) rmi4main.f90
rmi4main_FCFLAGS = -I../$(SRCOMMDIR)
rmi4main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librmi.la
librmi_la_SOURCES = $(SRC)
librmi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librmi_la_includedir = $(includedir)/$(SRINCDIR)
librmi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rmi.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)



