AC_INIT([slowroll], [0.1], [ringeval@gmail.com])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_RANLIB

FCFLAGS="-O -fopenmp"

SRCOMMDIR="common"
AC_SUBST([SRCOMMDIR])

FC_MODEXT="mod"
AC_SUBST([FC_MODEXT]) 

AM_INIT_AUTOMAKE(slowroll, 0.1)

AC_CONFIG_FILES([Makefile                 
		 src/common/Makefile
                 src/smallf/Makefile])
AC_OUTPUT
