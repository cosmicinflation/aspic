AC_INIT([slowroll], [0.1], [ringeval@gmail.com])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_RANLIB

FCFLAGS="-O -fopenmp"

SRCOMMDIR="common"
AC_SUBST([SRCOMMDIR])


AM_INIT_AUTOMAKE(slowroll, 0.1)

AC_CONFIG_FILES([Makefile                 
		 common/Makefile
                 smallf/Makefile])
AC_OUTPUT
