AC_INIT([aspic], [0.0.8], [chris@ringeval.com])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_RANLIB
AC_PROG_LIBTOOL

version_info="0:0:0"
AC_SUBST(version_info)

AC_LANG_PREPROC
FCFLAGS="-O"

SRCOMMDIR="common"
AC_SUBST([SRCOMMDIR])

SRINCDIR="aspic"
AC_SUBST([SRINCDIR])

FC_MODEXT="mod"
AC_SUBST([FC_MODEXT]) 

DATEXT="dat"
AC_SUBST([DATEXT]) 

AM_INIT_AUTOMAKE

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 man/Makefile
                 src/common/Makefile
                 src/lfi/Makefile
                 src/gmlfi/Makefile
                 src/rcmi/Makefile
                 src/rcqi/Makefile
                 src/pni/Makefile
                 src/mni/Makefile
                 src/esi/Makefile
                 src/pli/Makefile
                 src/kmii/Makefile
                 src/hf1i/Makefile
                 src/sfi/Makefile
                 src/li/Makefile
                 src/ii/Makefile
                 src/cwi/Makefile
                 src/kmiii/Makefile
                 src/hi/Makefile
                 src/twi/Makefile
                 src/mhi/Makefile
                 src/rpi/Makefile
                 src/dwi/Makefile
                 src/mssmi/Makefile
                 src/bsusybi/Makefile
                 src/ti/Makefile
                 src/bei/Makefile
                 src/psni/Makefile
                 src/ncki/Makefile
		 src/lmi/Makefile
                 src/rgi/Makefile
                 src/vhi/Makefile
                 src/dsi/Makefile
                 src/ssi/Makefile
                 src/oi/Makefile
                 src/ai/Makefile
                 src/cnai/Makefile
                 src/cnbi/Makefile
                 src/ripi/Makefile
                ])
AC_OUTPUT
