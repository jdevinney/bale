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

pkconfigdir=@libdir@/pkgconfig

include_HEADERS=std_options.h
lib_LIBRARIES=libstd_options.a

pkconfig_DATA = @PACKAGE_NAME@.pc
libstd_options_a_CPPFLAGS=@libdeps_CFLAGS@

common_sources=std_options.upc

if BUILD_UPC
libstd_options_a_SOURCES = $(common_sources:.c=__upc.upc)
endif

if BUILD_SHMEM
#libstd_options_a_LINK = $(LINK)
libstd_options_a_SOURCES = $(common_sources:.upc=__c.c)
libstd_options_a_CFLAGS = -DUSE_SHMEM
endif

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

