AC_INIT([srmodels], [0.1], [ringeval@iap.fr])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_RANLIB
AC_PROG_LIBTOOL

FCFLAGS="-O -fopenmp"

SRCOMMDIR="common"
AC_SUBST([SRCOMMDIR])

SRINCDIR="srmodels"
AC_SUBST([SRINCDIR])

FC_MODEXT="mod"
AC_SUBST([FC_MODEXT]) 

AM_INIT_AUTOMAKE(srmodels, 0.1)

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/common/Makefile
                 src/sfi/Makefile
                 src/lfi/Makefile
                 src/sfbi/Makefile
                 src/bi/Makefile
                 src/mlfi/Makefile
                 src/rmi/Makefile
                 src/vhi/Makefile
                ])
AC_OUTPUT
