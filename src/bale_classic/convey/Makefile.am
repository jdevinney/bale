# Copyright (c) 2020, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
#
# All rights reserved.
#
# This file is part of the conveyor package. For license information,
# see the LICENSE file in the top level directory of the distribution.


## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = $(warning_flags)

## Things that will be built and/or installed:

lib_LIBRARIES = libconvey.a
include_HEADERS = convey.h convey_alc8r.h convey_codec.h biconvey.h
if !WITH_UPC
  check_PROGRAMS = alltoallv histogram indexgather bigather meld
  convey_tuner = alltoallv
else
  # Because these programs will be linked differently, with $(UPCLINK),
  # we need to give them different names
  check_PROGRAMS = alltoallu histogupc indexgupc
  convey_tuner = alltoallu
endif
check_PROGRAMS += wringer
if WITH_SHMEM
  check_PROGRAMS += getput
endif
if WITH_SHMEM_PTR
  check_PROGRAMS += shmemptr
endif
dist_check_SCRIPTS = test_simple test_twohop test_tensor test_elastic \
    test_compress test_examples test_shmem_ptr

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = convey.pc

## Intermediate files produced by cases.pl

codelet_sources = accel_x.c
EXTRA_DIST = accel.c cases.pl
CLEANFILES = $(codelet_sources)

PERL = @PERL@
CASES = @CASES@
%_x.c: %.c
	$(PERL) $(CASES) -E $(srcdir) -p $< > $@


## The pieces of the library:

convey_headers = bolite.h convey.h convey_alc8r.h convey_codec.h convey_impl.h \
    common.h pivot.h porter.h porter_impl.h private.h router.h simple.h sorter.h tensor.h
convey_sources = circle.c codata.c common.c convey.c elastic.c packer.c porter.c \
    simple.c squeeze.c tensor.c trivial.c twohop.c $(codelet_sources)
convey_headers += biconvey.h biconvey_impl.h
convey_sources += biconvey.c bisimple.c bitensor.c

if WITH_SHMEM
  convey_sources += putport.c
endif
if WITH_MPI
  convey_sources += mpiport.c
endif
if WITH_UPC
  convey_sources += putport.c
endif

if !WITH_MPP_UTIL
if WITH_MPI
  convey_headers += mpp2mpi.h
  convey_sources += mpp2mpi.c
endif
if WITH_SHMEM
  convey_headers += mpp2shmem.h
  convey_sources += mpp2shmem.c
endif
if WITH_UPC
  convey_headers += mpp2upc.h
  convey_sources += mpp2upc.c
endif
if !WITH_MPI
if !WITH_SHMEM
if !WITH_UPC
  convey_headers += mpp2nil.h
  convey_sources += mpp2nil.c
endif
endif
endif
endif

libconvey_a_SOURCES = $(convey_headers)
if WITH_UPC
  nodist_libconvey_a_SOURCES = $(convey_sources:.c=__upc.upc)
  CLEANFILES += *__upc.upc
else
  libconvey_a_SOURCES += $(convey_sources)
endif

# Automatic linking to create UPC copies of C files.
%__upc.upc: %.c
	@echo "INFO: Linking UPC copy of   C file: $@"
	@rm -f $@
	@test -d $(dir $@) || $(MKDIR_P) $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; \
  then $(LN_S) $(notdir $<) $@ || cp -f $< $@; \
  else $(LN_S) $< $@ || cp -f $< $@; fi

## Benchmark and example programs:

if !WITH_UPC
  alltoallv_CFLAGS = $(AM_CFLAGS) $(ignore_deprecation)
  alltoallv_SOURCES = alltoallv.c a2a_align.c a2a_apps.c a2a_basic.c \
    alltoallv.h bolite.h common.h convey.h
  histogram_SOURCES = examples/histo.c examples/example.h bolite.h convey.h
  indexgather_SOURCES = examples/gather.c examples/example.h bolite.h convey.h
  bigather_SOURCES = examples/bigather.c examples/example.h bolite.h convey.h biconvey.h
  meld_SOURCES = examples/meld.c examples/example.h bolite.h convey.h biconvey.h
  alltoallv_LDADD = libconvey.a
  histogram_LDADD = libconvey.a
  indexgather_LDADD = libconvey.a
  bigather_LDADD = libconvey.a
  meld_LDADD = libconvey.a
else
  alltoallu_SOURCES = alltoallv.h bolite.h common.h convey.h
  histogupc_SOURCES = examples/example.h bolite.h convey.h
  indexgupc_SOURCES = examples/example.h bolite.h convey.h
  nodist_alltoallu_SOURCES = alltoallv__upc.upc \
    a2a_align__upc.upc a2a_apps__upc.upc a2a_basic__upc.upc
  nodist_histogupc_SOURCES = examples/histo__upc.upc
  nodist_indexgupc_SOURCES = examples/gather__upc.upc 
  alltoallu_LDADD = libconvey.a
  histogupc_LDADD = libconvey.a
  indexgupc_LDADD = libconvey.a
endif

## Testing:

wringer_SOURCES = wringer.c
wringer_LDADD = libconvey.a
if WITH_SHMEM
  getput_SOURCES = getput.c
endif
if WITH_SHMEM_PTR
  shmemptr_SOURCES = shmemptr.c
endif

LAUNCH ?= $(srcdir)/launchinfo
TESTS_ENVIRONMENT = LAUNCH="$(LAUNCH)"
TESTS =
if WITH_SHMEM_PTR
  TESTS += test_shmem_ptr
endif
TESTS += test_simple test_tensor test_elastic test_compress test_examples
if !WITH_UPC
  TESTS += test_twohop
endif

## Tuning:

tune: $(convey_tuner)
	@echo "Please wait a few minutes... (watch tune.log if you're bored)"
	LAUNCH="$(LAUNCH)" $(srcdir)/tune_tensor > tune.log
	$(AWK) '/args =/{b=($$6=="-B"); c=substr($$(b?7:6),3)} /bandw/{print b,c,$$NF*1e-6}' \
          tune.log | sort -k1n -k2n > tune.dat
	$(PERL) $(srcdir)/tune.pl < tune.dat

.PHONY: tune

## Distribution:

EXTRA_DIST += INSTALL LICENSE MAINTAIN README.md \
    launchinfo bench_simple.awk bench_tensor.awk tune_tensor tune.pl \
    narrative.h sorter.c
CLEANFILES += tests.out tune.log tune.dat

## Generate documentation using doxygen, if possible:

include m4/doxygen.am
