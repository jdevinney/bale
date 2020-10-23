#/******************************************************************
#
#
#  Copyright(C) 2020, Institute for Defense Analyses
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

bin_PROGRAMS = histo ig topo randperm permute_matrix transpose_matrix triangles sparse_matrix_io  sssp #unionfind

#BDEBUG = -g -O0
BDEBUG =

###################################################################
noinst_HEADERS = histo_src/histo.h 

histo_sources = histo_src/histo.upc 
histo_sources += histo_src/histo_agp.upc 
histo_sources += histo_src/histo_exstack.upc 
histo_sources += histo_src/histo_conveyor.upc 
histo_sources += histo_src/histo_exstack2.upc
histo_sources += histo_src/alternates/histo_exstack2_cyclic.upc
histo_sources += histo_src/alternates/histo_exstack2_goto.upc
histo_sources += histo_src/alternates/histo_exstack2_function.upc
histo_sources += histo_src/alternates/histo_exstack_function.upc

#@(info $$histo_sources is [${histo_sources}])

histo_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/histo_src 
histo_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
histo_SOURCES = $(histo_sources)
histo_LINK=$(UPCLINK)
endif

if BUILD_SHMEM
histo_SOURCES = $(histo_sources:.upc=__c.c) 
histo_CPPFLAGS += -DUSE_SHMEM
histo_LINK = $(LINK) $(BDEBUG)
endif

###################################################################
ig_sources = ig_src/ig.h 
ig_sources += ig_src/ig.upc 
ig_sources += ig_src/ig_agp.upc 
ig_sources += ig_src/ig_exstack.upc 
ig_sources += ig_src/ig_exstack2.upc 
ig_sources += ig_src/ig_conveyor.upc
ig_sources += ig_src/alternates/ig_exstack2_cyclic.upc
ig_sources += ig_src/alternates/ig_exstack2_goto.upc
ig_sources += ig_src/alternates/ig_exstack_function.upc
ig_sources += ig_src/alternates/ig_exstack_pkg.upc

ig_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/ig_src
ig_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
ig_SOURCES = $(ig_sources)
ig_LINK=$(UPCLINK)
endif

if BUILD_SHMEM
ig_SOURCES = $(ig_sources:.upc=__c.c)
ig_CPPFLAGS += -DUSE_SHMEM
ig_LINK=$(LINK) $(BDEBUG)
endif
###################################################################
topo_sources = topo_src/toposort.h
topo_sources += topo_src/toposort.upc
topo_sources += topo_src/toposort_agp.upc
topo_sources += topo_src/toposort_exstack.upc
topo_sources += topo_src/toposort_exstack2.upc
topo_sources += topo_src/toposort_conveyor.upc
topo_sources += topo_src/alternates/toposort_agp_oo.upc
topo_sources += topo_src/alternates/toposort_cooler.upc
topo_sources += topo_src/alternates/toposort_exstack_orig.upc

topo_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/topo_src
topo_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
topo_SOURCES = $(topo_sources)
topo_LINK=$(UPCLINK)
endif

if BUILD_SHMEM
topo_SOURCES = $(topo_sources:.upc=__c.c) 
topo_CPPFLAGS += -DUSE_SHMEM
topo_LINK=$(LINK) $(BDEBUG)
endif

###################################################################
noinst_HEADERS += randperm_src/randperm_alternates.h
randperm_sources = randperm_src/randperm.upc
randperm_sources += randperm_src/alternates/randperm_agp_opt.upc

randperm_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/randperm_src
randperm_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
randperm_SOURCES = $(randperm_sources)
randperm_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
randperm_SOURCES = $(randperm_sources:.upc=__c.c)
randperm_CPPFLAGS += -DUSE_SHMEM
randperm_LINK=$(LINK) $(BDEBUG)
endif

###################################################################
permute_matrix_sources = permute_matrix_src/permute_matrix.upc
permute_matrix_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/permute_matrix_src
permute_matrix_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
permute_matrix_SOURCES = $(permute_matrix_sources)
permute_matrix_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
permute_matrix_SOURCES = $(permute_matrix_sources:.upc=__c.c)
permute_matrix_CPPFLAGS += -DUSE_SHMEM
permute_matrix_LINK=$(LINK) $(BDEBUG)
endif

###################################################################
transpose_matrix_sources = transpose_matrix_src/transpose_matrix.upc
transpose_matrix_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/transpose_matrix_src
transpose_matrix_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
transpose_matrix_SOURCES = $(transpose_matrix_sources)
transpose_matrix_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
transpose_matrix_SOURCES = $(transpose_matrix_sources:.upc=__c.c)
transpose_matrix_CPPFLAGS += -DUSE_SHMEM
transpose_matrix_LINK=$(LINK) $(BDEBUG)
endif

###################################################################
noinst_HEADERS += triangle_src/triangle.h
triangles_sources = triangle_src/triangle.upc
triangles_sources += triangle_src/triangle_agp.upc
triangles_sources += triangle_src/triangle_exstack.upc
triangles_sources += triangle_src/triangle_exstack2.upc
triangles_sources += triangle_src/triangle_conveyor.upc
triangles_sources += triangle_src/alternates/triangle_agp_opt1.upc
triangles_sources += triangle_src/alternates/triangle_agp_opt2.upc
triangles_sources += triangle_src/alternates/triangle_agp_iter.upc
#triangles_sources += triangle_src/alternates/triangle_agp_oo.upc

triangles_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/triangle_src
triangles_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
triangles_SOURCES = $(triangles_sources)
triangles_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
triangles_SOURCES = $(triangles_sources:.upc=__c.c)
triangles_CPPFLAGS += -DUSE_SHMEM
triangles_LINK=$(LINK) $(BDEBUG)
endif

###################################################################
sparse_matrix_io_sources = sparse_matrix_io_src/sparse_matrix_io.upc

sparse_matrix_io_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/sparse_matrix_io_src
sparse_matrix_io_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
sparse_matrix_io_SOURCES = $(sparse_matrix_io_sources)
sparse_matrix_io_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
sparse_matrix_io_SOURCES = $(sparse_matrix_io_sources:.upc=__c.c)
sparse_matrix_io_CPPFLAGS += -DUSE_SHMEM
sparse_matrix_io_LINK=$(LINK) $(BDEBUG)
endif

###################################################################

###################################################################
noinst_HEADERS += sssp_src/sssp.h sssp_src/sssp_delta_common.h
sssp_sources = sssp_src/sssp.upc
sssp_sources += sssp_src/sssp_bellman_exstack.upc
sssp_sources += sssp_src/sssp_bellman_exstack2.upc
sssp_sources += sssp_src/sssp_bellman_conveyor.upc
sssp_sources += sssp_src/sssp_delta_exstack.upc
sssp_sources += sssp_src/sssp_delta_conveyor.upc
sssp_sources += sssp_src/sssp_delta_exstack2.upc
sssp_sources += sssp_src/sssp_delta_common.upc
sssp_sources += sssp_src/alternates/sssp_bellman_agp.upc
sssp_CPPFLAGS = $(BDEBUG) @libdeps_CFLAGS@ -I$(top_srcdir)/sssp_src
sssp_LDADD = @libdeps_LIBS@ -lm

if BUILD_UPC
sssp_SOURCES = $(sssp_sources)
sssp_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
sssp_SOURCES = $(sssp_sources:.upc=__c.c)
sssp_CPPFLAGS += -DUSE_SHMEM
sssp_LINK=$(LINK) $(BDEBUG)
endif

###################################################################
# Automatic linking as needed to create C copies of UPC files.               #
%__c.c: %.upc
	@echo $<
	@echo $@
	@pwd
	@echo "INFO: Linking   C copy of UPC file: $(realpath $<) to $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $(realpath $<) $@ || cp -f $< $@; fi


# Automatic linking as needed to create UPC copies of C files.               #
%__upc.upc: %.c
	@echo "INFO: Linking UPC copy of   C file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi
