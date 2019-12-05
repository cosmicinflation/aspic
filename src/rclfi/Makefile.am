SRC = rclficommon.f90 rclfi1sr.f90 rclfi1reheat.f90 \
	rclfi2sr.f90 rclfi2reheat.f90 \
	rclfi3sr.f90 rclfi3reheat.f90 \
	rclfi4sr.f90 rclfi4reheat.f90
MOD = rclficommon.$(FC_MODEXT)  \
	rclfi1sr.$(FC_MODEXT) rclfi1reheat.$(FC_MODEXT) \
	rclfi2sr.$(FC_MODEXT) rclfi2reheat.$(FC_MODEXT) \
	rclfi3sr.$(FC_MODEXT) rclfi3reheat.$(FC_MODEXT) \
	rclfi4sr.$(FC_MODEXT) rclfi4reheat.$(FC_MODEXT) 

check_PROGRAMS = rclfi1main rclfi2main rclfi3main rclfi4main
rclfi1main_SOURCES = $(SRC) rclfi1main.f90
rclfi1main_FCFLAGS = -I../$(SRCOMMDIR)
rclfi1main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
rclfi2main_SOURCES = $(SRC) rclfi2main.f90
rclfi2main_FCFLAGS = -I../$(SRCOMMDIR)
rclfi2main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
rclfi3main_SOURCES = $(SRC) rclfi3main.f90
rclfi3main_FCFLAGS = -I../$(SRCOMMDIR)
rclfi3main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a
rclfi4main_SOURCES = $(SRC) rclfi4main.f90
rclfi4main_FCFLAGS = -I../$(SRCOMMDIR)
rclfi4main_LDADD = ../$(SRCOMMDIR)/libsrcommon.a

noinst_LTLIBRARIES = librclfi.la
librclfi_la_SOURCES = $(SRC)
librclfi_la_FCFLAGS = -I../$(SRCOMMDIR) $(AM_FCFLAGS)
librclfi_la_includedir = $(includedir)/$(SRINCDIR)
librclfi_la_include_HEADERS = $(MOD)

man_MANS = aspic_rclfi1.3 aspic_rclfi2.3 aspic_rclfi3.3 aspic_rclfi4.3

clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)


.NOTPARALLEL:
