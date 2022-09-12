AC_INIT([aspic],[1.2.0],[aspic@ringeval.info],[],[http://curl.irmp.ucl.ac.be/~chris/aspic.html])

AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_FC
AC_PROG_CPP
LT_INIT


version_info="4:0:0"
AC_SUBST(version_info)

AC_LANG_PREPROC
FCFLAGS+=" "

SRCOMMDIR="common"
AC_SUBST([SRCOMMDIR])

AC_ARG_ENABLE(quad-precision,
	[AS_HELP_STRING([--enable-quad-precision],[compile aspic in quadruple precision])],
	ok=$enableval, ok=no)
AM_CONDITIONAL([USE_QUAD],[test "$ok" = "yes"])
AM_COND_IF([USE_QUAD],
[
	FCFLAGS+=" -DQUADPREC"
	PREC_SUFFIX=q
	SRINCDIR="aspicq"
	AC_MSG_WARN([attempting to compile in quadruple precision!])
],
[
	SRINCDIR="aspic"
]
)


AC_ARG_ENABLE(openmp,
	[AS_HELP_STRING([--enable-openmp],[enable openmp parallel instructions in aspic])],
	ok=$enableval, ok=yes)
AM_CONDITIONAL([USE_OPENMP],[test "$ok" = "yes"])
AM_COND_IF([USE_OPENMP],
[
	FCFLAGS+=" -fopenmp"
],
[
	AC_MSG_WARN([disabling openmp parallel instructions!])
	FCFLAGS+=" -DNOOMP"
]
)

AC_SUBST(PREC_SUFFIX)
AC_SUBST([SRINCDIR])

FC_MODEXT="mod"
AC_SUBST([FC_MODEXT]) 

DATEXT="dat"
AC_SUBST([DATEXT]) 

AC_CONFIG_TESTDIR([tests])

AC_CONFIG_FILES([Makefile
		 tests/Makefile
                 src/Makefile
                 man/Makefile
		 src/common/Makefile
                 src/lfi/Makefile
                 src/gmlfi/Makefile
                 src/rcmi/Makefile
                 src/rcqi/Makefile
                 src/ni/Makefile
                 src/esi/Makefile
                 src/pli/Makefile
                 src/kmii/Makefile
                 src/hf1i/Makefile
                 src/sfi/Makefile
                 src/li/Makefile
                 src/ii/Makefile
                 src/cwi/Makefile
                 src/kmiii/Makefile
                 src/si/Makefile
		 src/hi/Makefile
                 src/twi/Makefile
                 src/mhi/Makefile
                 src/rpi/Makefile
                 src/dwi/Makefile
                 src/gmssmi/Makefile
                 src/bsusybi/Makefile
                 src/ti/Makefile
                 src/bei/Makefile
                 src/psni/Makefile
                 src/ncki/Makefile
		 src/lmi/Makefile
                 src/rgi/Makefile
                 src/vhi/Makefile
                 src/dsi/Makefile
                 src/ssbi/Makefile
                 src/oi/Makefile
                 src/ai/Makefile
                 src/cnai/Makefile
                 src/cnbi/Makefile
                 src/oripi/Makefile
                 src/csi/Makefile
                 src/sbi/Makefile
                 src/cnci/Makefile
                 src/rmi/Makefile
                 src/lpi/Makefile
                 src/cndi/Makefile
                 src/rchi/Makefile
                 src/osti/Makefile
                 src/wri/Makefile
                 src/imi/Makefile
                 src/gripi/Makefile
                 src/bi/Makefile
                 src/kklti/Makefile
                 src/nfi/Makefile
                 src/ccsi/Makefile
		 src/di/Makefile
		 src/hni/Makefile
		 src/ncli/Makefile
		 src/vfmi/Makefile
		 src/fi/Makefile
		 src/ahi/Makefile
		 src/sbki/Makefile
		 src/sdi/Makefile
		 src/saai/Makefile
		 src/sabi/Makefile
		 src/saci/Makefile
		 src/hbi/Makefile
		 src/shi/Makefile
		 src/rcipi/Makefile
 		 src/gdwi/Makefile
		 src/saii/Makefile
		 src/saiii/Makefile
		 src/rclfi/Makefile
		 src/dei/Makefile
		 src/pai/Makefile
                ])
AC_OUTPUT
