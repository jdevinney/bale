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

common_sources = exstack.upc 
pkconfigdir=@libdir@/pkgconfig

include_HEADERS=exstack.h
lib_LIBRARIES = libexstack.a

pkconfig_DATA = @PACKAGE_NAME@.pc

libexstack_a_CPPFLAGS=

# handle UPC library
if BUILD_UPC
libexstack_a_SOURCES = $(common_sources:.c=__upc.upc)
libexstack_a_SOURCES += exstack2.upc
endif

if BUILD_SHMEM
libexstack_a_SOURCES = $(common_sources:.upc=__c.c)
libexstack_a_SOURCES += exstack2_shmem.c
libexstack_a_CFLAGS = -DUSE_SHMEM 
endif

libexstack_a_CPPFLAGS += @CPPFLAGS@ @libdeps_CFLAGS@

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


