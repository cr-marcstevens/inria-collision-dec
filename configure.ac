AC_INIT([ISD], [1.0], [collision-dec-bugs@lists.gforge.inria.fr])
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects]) #subdir-objects creates objects in sub-directories
AC_PROG_CC
AC_PROG_CC_STDC #to have a C99 compiler
AC_CONFIG_HEADERS([config.h])
AM_PROG_CC_C_O
AC_PROG_RANLIB

AM_PATH_PYTHON([2.7])

PKG_CHECK_MODULES([M4RI], [m4ri >= 20130416])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
