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

common_sources = libgetput.upc
pkconfigdir=@libdir@/pkgconfig
noinst_HEADERS=knuth_rng_double_2019.h
include_HEADERS=libgetput.h
#lib_LIBRARIES = 
lib_LIBRARIES = liblibgetput.a

pkconfig_DATA = @PACKAGE_NAME@.pc

if BUILD_UPC
liblibgetput_a_SOURCES = $(common_sources:.c=__upc.upc)
endif

if BUILD_SHMEM
#liblibgetput_a_LINK = $(LINK)
liblibgetput_a_SOURCES = $(common_sources:.upc=__c.c)
liblibgetput_a_CFLAGS = -DUSE_SHMEM -D_XOPEN_SOURCE=700
endif

# Automatic linking as needed to create C copies of UPC files.               #
%__c.c: %.upc
	@echo "INFO: Linking   C copy of UPC file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi


