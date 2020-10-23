#/******************************************************************
#
#
#  Copyright(C) 2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
#
#  All rights reserved.
#  
#  This file is a part of Bale.  For license information see the
#  LICENSE file in the top level directory of the distribution.
#  
# *****************************************************************/ 

common_sources  = spmat_agp.upc 
common_sources += spmat_exstack.upc
common_sources += spmat_exstack2.upc 
common_sources += spmat_conveyor.upc
common_sources += spmat_io.upc
common_sources += spmat_utils.upc
common_sources += geometric.upc
pkconfigdir=@libdir@/pkgconfig

include_HEADERS=spmat.h spmat_enums.h

lib_LIBRARIES = libspmat.a

pkconfig_DATA = @PACKAGE_NAME@.pc
libspmat_a_CPPFLAGS = -g -O2
#libspmat_a_CPPFLAGS = -g -O0

if BUILD_UPC
libspmat_a_SOURCES = $(common_sources:.c=__upc.upc)
endif

if BUILD_SHMEM
#libspmat_a_LINK = $(LINK)
libspmat_a_SOURCES = $(common_sources:.upc=__c.c)
libspmat_a_CPPFLAGS += -DUSE_SHMEM
endif

libspmat_a_CPPFLAGS += @CPPFLAGS@ @libdeps_CFLAGS@

# Automatic linking as needed to create C copies of UPC files.               #
%__c.c: %.upc
	@echo "INFO: Linking   C copy of UPC file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi


# Automatic linking as needed to create UPC copies of C files.               #
%__upc.upc: %.c
	@echo "INFO: Linking UPC copy of   C file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi


