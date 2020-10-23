#/******************************************************************
#
#
#  Copyright(C) 2019-2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
#
#  All rights reserved.
#  
#  This file is a part of Bale.  For license information see the
#  LICENSE file in the top level directory of the distribution.
#  
# *****************************************************************/ 

AC_PREREQ([2.60])     #Require Minimum Autotools version 2.60

AC_INIT([spmat],   [2.1.0],   [bale@super.org])
AM_INIT_AUTOMAKE([1.10 foreign no-define tar-ustar])
##LT_INIT ## enable libtool

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([../])

AC_ARG_WITH([upc],
            AS_HELP_STRING([--with-upc], [build a UPC version]),
            [], [with_upc=yes])
AC_ARG_WITH([shmem],
            AS_HELP_STRING([--with-shmem], [build a SHMEM version]),
            [], [with_shmem=no])

AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CC_C99

mpp_use_shmem=no
if test "x$with_shmem" = xyes; then
   with_upc=no
   mpp_use_shmem=yes
fi

if test "x$with_upc" = xyes; then
AM_PROG_UPC
else
  # must define this ugly conditional to keep automake happy
  AM_CONDITIONAL(am__fastdepUPC, [test])
fi

AS_IF([test "x$mpp_use_shmem" = xyes],
[AC_SEARCH_LIBS([shmemx_alltoallv], [],
    [AC_DEFINE([HAVE_SHMEMX_ALLTOALLV], [1], [Define if shmem_alltoallv exists.])])
AC_SEARCH_LIBS([shmemx_team_alltoallv], [],
    [AC_DEFINE([HAVE_SHMEMX_TEAM_ALLTOALLV], [1], [Define if shmem_team_alltoallv exists.])])
AC_SEARCH_LIBS([shmemx_putmem_signal], [],
    [AC_DEFINE([HAVE_SHMEMX_PUTMEM_SIGNAL], [1], [Define if shmem_putmem_signal exists.])])
AC_SEARCH_LIBS([shmem_free], [],
    [AC_DEFINE([HAVE_SHMEM_FREE], [1], [Define if shmem_free exists.])])
])

AM_CONDITIONAL(BUILD_UPC, [test "x$with_upc" = xyes])
AM_CONDITIONAL(BUILD_SHMEM, [test "x$mpp_use_shmem" = xyes])

PKG_CHECK_MODULES([libdeps],
  [
        exstack                
        convey
        libgetput
  ]
)

AC_CONFIG_FILES([
spmat.pc:spmat.pc.in
Makefile
])

AC_OUTPUT

