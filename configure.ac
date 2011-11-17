AC_INIT([srmodels], [0.0.3], [ringeval@iap.fr])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_RANLIB
AC_PROG_LIBTOOL

version_info="0:0:0"
AC_SUBST(version_info)

FCFLAGS="-O -fopenmp"

SRCOMMDIR="common"
AC_SUBST([SRCOMMDIR])

SRINCDIR="srmodels"
AC_SUBST([SRINCDIR])

FC_MODEXT="mod"
AC_SUBST([FC_MODEXT]) 

DATEXT="dat"
AC_SUBST([DATEXT]) 

AM_INIT_AUTOMAKE

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/common/Makefile
                 src/lfi/Makefile
                 src/mlfi/Makefile
                 src/rcmi/Makefile
                 src/rcqi/Makefile
                 src/pni/Makefile
                 src/mni/Makefile
                 src/esi/Makefile
                 src/pli/Makefile
                 src/kmii/Makefile
                 src/hf1i/Makefile
                ])
AC_OUTPUT
