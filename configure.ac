AC_INIT([ISD], [1.0], [collision-dec-bugs@lists.gforge.inria.fr])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC
AC_PROG_CC_STDC #to have a C99 compiler
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects]) #create objects in sub-directories
AM_PROG_CC_C_O
AC_PROG_RANLIB


# Get m4ri library and include locations
AC_ARG_WITH([m4ri-include-path],
  [AS_HELP_STRING([--with-m4ri-include-path],
    [location of the m4ri headers, defaults to /usr/include/m4ri])],
  [M4RI_INCLUDE="$withval"],
  [M4RI_INCLUDE="/usr/include/m4ri"])
AC_SUBST([M4RI_INCLUDE])

AC_ARG_WITH([m4ri-lib-path],
  [AS_HELP_STRING([--with-m4ri-lib-path],
		[location of the m4ri librarie, default to /usr/lib])],
  [M4RI_LIB="$withval"],
  [M4RI_LIB="/usr/lib"])
AC_SUBST([M4RI_LIB])


AC_CONFIG_FILES([Makefile src/Makefile src/cpucycles/Makefile])
AC_OUTPUT
