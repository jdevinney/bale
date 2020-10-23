# Copyright (c) 2020, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
#
# All rights reserved.
#
# This file is part of the conveyor package. For license information,
# see the LICENSE file in the top level directory of the distribution.


## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
AC_INIT([convey], [0.6.0], [bale@super.org])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([simple.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign subdir-objects])

# Configure options and variables
AC_ARG_VAR([PERL], [Perl interpreter])
AC_ARG_VAR([CASES], [Perl script for preprocessing source code])
AC_ARG_VAR([DOXYGEN], [Doxygen documentation generation program])
AC_ARG_VAR([MSGBYTES], [Preferred message size in bytes])
AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi], [build on MPI alone])],
            [], [with_mpi=no])
AC_ARG_WITH([mpp],
            [AS_HELP_STRING([--with-mpp], [build atop mpp_utilV4])],
            [], [with_mpp=no])
AC_ARG_WITH([shmem],
            [AS_HELP_STRING([--with-shmem], [build atop SHMEM])],
            [], [with_shmem=no])
AC_ARG_WITH([upc],
            [AS_HELP_STRING([--with-upc], [build a UPC version])],
            [], [with_upc=no])
AC_ARG_WITH([shmem-ptr],
            [AS_HELP_STRING([--without-shmem-ptr], [avoid use of shmem_ptr()])],
            [], [with_shmem_ptr=yes])
AC_ARG_ENABLE([profiling],
    [AS_HELP_STRING([--enable-profiling], [enable mpp_util profiling])],
    [], [enable_profiling=no])
AC_ARG_ENABLE([nonblocking],
    [AS_HELP_STRING([--disable-nonblocking], [disable nonblocking puts])],
    [], [enable_nonblocking=yes])

if   test "x$with_shmem" = xyes -a "x$with_mpi" = xyes; then
    AC_MSG_ERROR([--with-shmem and --with-mpi are incompatible])
elif test "x$with_shmem" = xyes -a "x$with_mpp" = xyes; then
    AC_MSG_ERROR([--with-shmem and --with-mpp are incompatible])
elif test "x$with_shmem" = xyes -a "x$with_upc" = xyes; then
    AC_MSG_ERROR([--with-shmem and --with-upc are incompatible])
elif test "x$with_mpp" = xyes -a "x$with_upc" = xyes; then
    AC_MSG_ERROR([--with-mpp and --with-upc are incompatible])
elif test "x$with_mpp" = xyes -a "x$with_mpi" = xyes; then
    AC_MSG_ERROR([--with-mpp and --with-mpi are incompatible])
elif test "x$with_upc" = xyes -a "x$with_mpi" = xyes; then
    AC_MSG_ERROR([--with-upc and --with-mpi are incompatible])
fi

if test "x$with_mpp" = xno -a "x$enable_profiling" = xyes; then
    AC_MSG_ERROR([--enable-profiling requires --with-mpp])
fi

REQUIREMENTS=

# Checks for programs.  (Before we can use AC_CHECK_HEADERS, we
# must use AC_PROG_CPP to find the C preprocessor.)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_PATH_PROG([PERL], [perl])
if ! test "$PERL" ; then
    AC_MSG_ERROR([perl not found; use PERL=... or update PATH])
fi
if ! test "$CASES" ; then
    CASES='$(srcdir)/cases.pl'
fi
AC_PATH_PROG([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
    AC_MSG_WARN([doxygen not found; continuing without it])
fi

warning_flags=
ignore_deprecation=
if test "x$with_upc" = xyes; then
  AM_PROG_UPC
  AC_PROG_MKDIR_P
  AC_PROG_LN_S
else
  # must define this ugly conditional to keep automake happy
  AM_CONDITIONAL(am__fastdepUPC, [test])
  AX_CHECK_COMPILER_FLAGS([-Wall -Wno-parentheses -Wno-unused -Wwrite-strings],
    [warning_flags="-Wall -Wno-parentheses -Wno-unused -Wwrite-strings"], [])
  AX_CHECK_COMPILER_FLAGS([-Wno-deprecated-declarations],
    [ignore_deprecation="-Wno-deprecated-declarations"], [])
fi

if test "x$with_mpp" = xyes; then
    ACX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS"
                 CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
                 CC="$PTHREAD_CC"],
                [AC_MSG_ERROR([cannot determine how to compile with pthreads])])
fi

# Define how to push and pop the UPC language
AC_DEFUN([AX_LANG_PUSH_UPC], [AC_LANG_PUSH([C])
ac_compile='$UPC -c $UPCFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$UPC -o conftest$ac_exeext $UPCFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=no])

AC_DEFUN([AX_LANG_POP_UPC], [AC_LANG_POP([C])])

# Consider finding mpp_util via pkg-config...
# PKG_CHECK_MODULES([MPP_UTIL], [mpp_utilV4], , )
# CPPFLAGS="$CPPFLAGS $MPP_UTIL_CFLAGS"
# LIBS="$MPP_UTIL_LIBS $LIBS"

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADERS([inttypes.h stdalign.h])
AC_TYPE_UINTPTR_T
if test "x$ac_cv_c_uintptr_t" = xno; then
    AC_MSG_ERROR([type uintptr_t is not available])
fi
if test "x$with_upc" = xno; then
    AC_CHECK_HEADERS([stdatomic.h])
    AC_CHECK_TYPES([_Atomic uint64_t], [], [], [[#include <stdatomic.h>
#include <stdint.h>]])
    AC_MSG_CHECKING([whether atomics are usable])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdatomic.h>]],
[[_Atomic long long foo; atomic_store(&foo, 1);]])],
        [AC_MSG_RESULT([yes])
         AC_DEFINE([HAVE_USABLE_ATOMICS], [1],
                   [Define to 1 if _Atomic long long and atomic_store() are usable.])],
        [AC_MSG_RESULT([no])])
fi

# Checks for header files
AC_CHECK_HEADERS([unistd.h], ,
    [AC_MSG_ERROR([Couldn't find unistd.h. This is bad.])])
if test "x$with_mpi" = xyes; then
    AC_CHECK_HEADERS([mpi.h], ,
        [AC_MSG_ERROR([Couldn't find mpi.h. Update CC or CPPFLAGS and try again.])])
fi
if test "x$with_shmem" = xyes; then
    AC_CHECK_HEADERS([shmem.h], ,
        [AC_MSG_ERROR([Couldn't find shmem.h. Update CPPFLAGS and try again.])])
fi
## Set up mpp_utilV4 in the conventional way for now
if test "x$with_mpp" = xyes; then
    AC_CHECK_HEADERS([micro64.h bitops.h], ,
        [AC_MSG_ERROR([Couldn't find bitops.h and micro64.h. Update CPPFLAGS and try again.])])
    AC_CHECK_HEADERS([actimer.h], ,
        [AC_MSG_ERROR([Couldn't find actimer.h. Update CPPFLAGS and try again.])],
        [#define ACTIMER_MODULE_NAME conftest
#define ACTIMER_MAIN])
    AC_CHECK_HEADERS([mpp_utilV4.h], ,
        [AC_MSG_ERROR([Couldn't find mpp_utilV4.h. Update CPPFLAGS and try again.])])
fi

# Checks for dependencies
AC_CHECK_FUNCS([gethostname])
AC_CHECK_LIB([m], [sqrt], ,
    [AC_MSG_ERROR([Couldn't find the math library. Something is very wrong.])])
if test "x$with_mpi" = xyes; then
    AC_CHECK_FUNC([MPI_Init], ,
        [AC_MSG_ERROR([Couldn't find MPI_Init. Update LIBS, LDFLAGS, or CC and try again.])])
fi
if test "x$with_shmem" = xyes; then
    AC_CHECK_FUNC([shmem_my_pe], ,
        [AC_MSG_ERROR([Couldn't find shmem_my_pe. Update LIBS, LDFLAGS, or CC and try again.])])
fi
if test "x$with_mpp" = xyes; then
    AC_CHECK_LIB([rt], [clock_gettime], ,
        [AC_MSG_ERROR([Couldn't find the real-time library. Something is very wrong.])])
    AC_CHECK_LIB([actimer], [actimer_summarize], ,
        [AC_MSG_ERROR([Couldn't find libactimer. Update LDFLAGS and try again.])])
fi

# Check behavior of getopt() [Linux vs. BSD]
AC_LANG_PUSH([C])
AC_MSG_CHECKING([for optreset in getopt()])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <unistd.h>]], [return optreset;])],
    [AC_MSG_RESULT([yes])
     have_optreset=yes],
    [AC_MSG_RESULT([no])
     have_optreset=no])
AC_LANG_POP([C])

mpp_use_mpi=no
mpp_use_shmem=no
if test "x$with_mpp" = xyes; then
    REQUIREMENTS="mpp_utilV4, actimer"
    if test "x$enable_profiling" = xyes; then
        AC_CHECK_LIB([z], [gzopen], ,
            [AC_MSG_ERROR([Couldn't find libz. Something is very wrong.])])
        AC_CHECK_LIB([mpp_utilV4_prof], [mpp_alltoall], ,
            [AC_MSG_ERROR([Couldn't find libmpp_utilV4_prof. Update LDFLAGS and try again.])])
        
    else
        AC_CHECK_LIB([mpp_utilV4], [mpp_alltoall], ,
            [AC_MSG_ERROR([Couldn't find libmpp_utilV4. Update LDFLAGS and try again.])])
    fi
    AC_LANG_PUSH([C])
    AC_MSG_CHECKING([for MPI mode])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpp_utilV4.h>
#if !MPP_USE_MPI
#error MPP_USE_MPI is false
#endif]], [return 0;])],
        [AC_MSG_RESULT([yes])
         mpp_use_mpi=yes],
        [AC_MSG_RESULT([no])])
    AC_MSG_CHECKING([for SHMEM mode])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpp_utilV4.h>
#if !MPP_USE_SHMEM
#error MPP_USE_SHMEM is false
#endif]], [return 0;])],
        [AC_MSG_RESULT([yes])
         mpp_use_shmem=yes],
        [AC_MSG_RESULT([no])])
    AC_LANG_POP([C])
elif test "x$with_upc" = xyes; then
    AX_LANG_PUSH_UPC
    AC_CHECK_HEADERS([upc.h])
    AC_CHECK_HEADERS([upc_castable.h], [],
        [AC_DEFINE([HAVE_UPC_CASTABLE_H], [1], [Define if upc_castable.h exists.])],
        [#include <upc.h>])
    AX_LANG_POP_UPC
fi

AS_IF([test "x$mpp_use_shmem" = xyes -o "x$with_shmem" = xyes],
[AC_CHECK_HEADERS([shmem.h])
AC_SEARCH_LIBS([shmemx_alltoallv], [],
    [AC_DEFINE([HAVE_SHMEMX_ALLTOALLV], [1], [Define if shmemx_alltoallv exists.])])
AC_SEARCH_LIBS([shmemx_team_alltoallv], [],
    [AC_DEFINE([HAVE_SHMEMX_TEAM_ALLTOALLV], [1], [Define if shmemx_team_alltoallv exists.])])
AC_SEARCH_LIBS([shmemx_putmem_signal], [],
    [AC_DEFINE([HAVE_SHMEMX_PUTMEM_SIGNAL], [1], [Define if shmemx_putmem_signal exists.])])
AC_SEARCH_LIBS([shmem_putmem_nbi], [],
    [AC_DEFINE([HAVE_SHMEM_PUTMEM_NBI], [1], [Define if shmem_putmem_nbi exists.])])
AC_SEARCH_LIBS([shmem_free], [],
    [AC_DEFINE([HAVE_SHMEM_FREE], [1], [Define if shmem_free exists.])])
AC_SEARCH_LIBS([shmem_align], [],
    [AC_DEFINE([HAVE_SHMEM_ALIGN], [1], [Define if shmem_align exists.])])
AC_SEARCH_LIBS([shmem_global_exit], [],
    [AC_DEFINE([HAVE_SHMEM_GLOBAL_EXIT], [1], [Define if shmem_global_exit exists.])])
], [with_shmem_ptr=no])
AS_IF([test "x$with_shmem_ptr" = xyes],
    [AC_SEARCH_LIBS([shmem_ptr], [], [], [with_shmem_ptr=no])])

AS_IF([test "x$with_mpp" = xyes],
    [AC_DEFINE([HAVE_MPP_UTIL], [1], [Define if the build is based on mpp_utilV4.])])
AS_IF([test "x$with_upc" = xyes],
    [AC_DEFINE([MPP_USE_UPC], [1], [Define if the build is based on UPC.])])
AS_IF([test "x$with_mpi" = xyes],
    [AC_DEFINE([MPP_RAW_MPI], [1], [Define if the build is based on raw MPI.])])
AS_IF([test "x$with_shmem" = xyes],
    [AC_DEFINE([MPP_RAW_SHMEM], [1], [Define if the build is based on raw SHMEM.])])
AS_IF([test "x$with_shmem_ptr" = xyes],
    [AC_DEFINE([HAVE_SHMEM_PTR], [1], [Define if shmem_ptr exists and should be used.])])
AS_IF([test "x$enable_profiling" = xyes],
    [AC_DEFINE([ENABLE_PROFILING], [1], [Define to activate mpp profiling and actimer.])])
AS_IF([test "x$have_optreset" = xyes],
    [AC_DEFINE([HAVE_OPTRESET], [1], [Define if getopt() understands optreset.])])

# Decide what codelets to build
AC_LANG_PUSH([C])
AC_MSG_CHECKING([for SSSE3 support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <tmmintrin.h>
#ifndef __SSSE3__
#error __SSSE3__ is not defined
#endif]],
        [__m128i x = _mm_setzero_si128(); x = _mm_alignr_epi8(x,x,1); return _mm_cvtsi128_si64(x);])],
    [AC_MSG_RESULT([yes])
     enable_ssse3=yes],
    [AC_MSG_RESULT([no])
     enable_ssse3=no])
AC_MSG_CHECKING([for AVX2 support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <immintrin.h>
#ifndef __AVX2__
#error __AVX2__ is not defined
#endif]],
        [__m256i x = _mm256_setzero_si256(); return _mm256_movemask_epi8(_mm256_add_epi32(x,x));])],
    [AC_MSG_RESULT([yes])
     enable_avx2=yes],
    [AC_MSG_RESULT([no])
     enable_avx2=no])
AC_MSG_CHECKING([for BMI2 support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <immintrin.h>
#ifndef __BMI2__
#error __BMI2__ is not defined
#endif]],
        [return _pext_u64(0,0);])],
    [AC_MSG_RESULT([yes])
     enable_bmi2=yes],
    [AC_MSG_RESULT([no])
     enable_bmi2=no])
AC_LANG_POP([C])
AS_IF([test "x$enable_ssse3" = xyes],
   [AC_DEFINE([HAVE_SSSE3], [1], [SSSE3 instructions available])])
AS_IF([test "x$enable_avx2" = xyes],
   [AC_DEFINE([HAVE_AVX2], [1], [AVX2 instructions available])])
AS_IF([test "x$enable_bmi2" = xyes],
   [AC_DEFINE([HAVE_BMI2], [1], [BMI2 instructions available])])

# Define conditions for automake
AM_CONDITIONAL(WITH_MPI, [test "x$mpp_use_mpi" = xyes -o "x$with_mpi" = xyes])
AM_CONDITIONAL(WITH_MPP_UTIL, [test "x$with_mpp" = xyes])
AM_CONDITIONAL(WITH_SHMEM, [test "x$mpp_use_shmem" = xyes -o "x$with_shmem" = xyes])
AM_CONDITIONAL(WITH_UPC, [test "x$with_upc" = xyes])
AM_CONDITIONAL(WITH_SHMEM_PTR, [test "x$with_shmem_ptr" = xyes])
AM_CONDITIONAL(HAVE_DOXYGEN, [test -n "$DOXYGEN"])

# Default parameters
if ! test "$MSGBYTES" ; then
    MSGBYTES=10000
fi
AC_DEFINE_UNQUOTED([CONVEY_BUFFER_SIZE], [$MSGBYTES],
    [Preferred message size, or buffer capacity, in bytes])
if test \( "x$with_shmem" = xno -a "x$mpp_use_shmem" = xno \) -o \
    "x$ac_cv_search_shmem_putmem_nbi" = xno ; then
    enable_nonblocking=no
fi
if test "x$enable_nonblocking" = xyes ; then
  enb_flag=1
else
  enb_flag=0
fi
AC_DEFINE_UNQUOTED([ENABLE_NONBLOCKING], [$enb_flag], [Define to enable nonblocking puts])

# Produce output
AC_SUBST([VERSION])
AC_SUBST([REQUIREMENTS])
AC_SUBST([warning_flags])
AC_SUBST([ignore_deprecation])
AC_CONFIG_FILES([Makefile convey.pc])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES(Doxyfile)])
AC_OUTPUT

# Report the results
if test "x$with_upc" = xyes; then
    prog_model=UPC
elif test "x$with_mpi" = xyes; then
    prog_model=MPI
elif test "x$with_shmem" = xyes; then
    prog_model=SHMEM
elif test "x$mpp_use_mpi" = xyes; then
    prog_model=MPP_UTIL/MPI
elif test "x$mpp_use_shmem" = xyes; then
    prog_model=MPP_UTIL/SHMEM
elif test "x$with_mpp" = xyes; then
    prog_model=MPP_UTIL/NONE
else
    prog_model=C
fi

echo "
     $PACKAGE_NAME $PACKAGE_VERSION
     Prefix.............: $prefix
     Programming Model..: $prog_model
     Message size.......: $MSGBYTES bytes
     Nonblocking puts...: $enable_nonblocking
     Using shmem_ptr....: $with_shmem_ptr
     MPP Profiling......: $enable_profiling
     SSSE3 Codelets.....: $enable_ssse3
     AVX2 Codelets......: $enable_avx2
     BMI2 Codelets......: $enable_bmi2
     Cases script.......: $PERL $CASES
     C Compiler.........: $CC $CFLAGS $CPPFLAGS
     Linker.............: $CC $LDFLAGS $LIBS
"
