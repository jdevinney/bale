#/******************************************************************
#
#
#  Copyright(C) 2019-2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# 
#
#  All rights reserved.
#  
#   This file is a part of Bale.  For license information see the
#   LICENSE file in the top level directory of the distribution.
#  
# 
# *****************************************************************/ 

AC_PREREQ([2.60])                   #Require Minimum Autotools version 2.60

AC_INIT([apps],   [2.1.0],   [bale@super.org])

AM_INIT_AUTOMAKE([1.10 foreign no-define tar-ustar subdir-objects])

AC_ARG_WITH([upc],
            AS_HELP_STRING([--with-upc], [build a UPC version]),
            [], [with_upc=yes])
AC_ARG_WITH([shmem],
            AS_HELP_STRING([--with-shmem], [build a SHMEM version]),
            [], [with_shmem=no])

# check for programs and compilers
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

AM_CONDITIONAL(BUILD_UPC, [test "x$with_upc" = xyes])
AM_CONDITIONAL(BUILD_SHMEM, [test "x$mpp_use_shmem" = xyes])

# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_COMPILE_FLAG(FLAG, [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   Check whether the given FLAG works with the current language's compiler
#   or gives an error.  (Warnings, however, are ignored)
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   INPUT gives an alternative input source to AC_COMPILE_IFELSE.
#
#   NOTE: Implementation based on AX_CFLAGS_GCC_OPTION. Please keep this
#   macro in sync with AX_CHECK_{PREPROC,LINK}_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_IF
AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

AX_CHECK_COMPILE_FLAG([-std=c11], [
  AS_VAR_APPEND([CFLAGS]," -std=c11")
],
[
AX_CHECK_COMPILE_FLAG([-hstd=c11],[  AS_VAR_APPEND([CFLAGS]," -hstd=c11")],
[echo "C compiler cannot compile C11 code"
exit -1])
])

PKG_CHECK_MODULES([libdeps],
  [
        std_options
        spmat
        convey
        exstack
        libgetput
  ]
)

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

# config files
AC_CONFIG_FILES([
Makefile
])

AC_OUTPUT
