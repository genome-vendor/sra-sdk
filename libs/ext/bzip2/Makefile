# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================


default: std

TOP ?= $(shell ../../../build/abspath.sh ../../..)
MODULE = libs/ext/bzip2

include $(TOP)/build/Makefile.env

INT_LIBS = \

EXT_LIBS = \
	libbz2

ALL_LIBS = \
	$(INT_LIBS) \
	$(EXT_LIBS)

#-------------------------------------------------------------------------------
# outer targets
#
std: makedirs download-$(OS)
	@ $(MAKE_CMD) $(TARGDIR)/$@

all: makedirs download
	@ $(MAKE_CMD) $(TARGDIR)/$@

$(INT_LIBS): makedirs download-$(OS)
	@ $(MAKE_CMD) $(ILIBDIR)/$@

$(EXT_LIBS): makedirs download-$(OS)
	@ $(MAKE_CMD) $(LIBDIR)/$@

.PHONY: all std $(ALL_LIBS)

#-------------------------------------------------------------------------------
# all
#
$(TARGDIR)/all: $(LIBDIR)/libbz2

.PHONY: $(TARGDIR)/all

#-------------------------------------------------------------------------------
# std
#
$(TARGDIR)/std: $(LIBDIR)/libbz2-$(OS)

.PHONY: $(TARGDIR)/std

#-------------------------------------------------------------------------------
# clean
#
clean: stdclean

deep-clean: clean
	rm -f *.[ch] *.def *.vers *.gz

.PHONY: clean deep-clean

#-------------------------------------------------------------------------------
# libbz2-$(OS)
#
$(LIBDIR)/libbz2-linux: $(LIBDIR)/libbz2-linux-$(BITS)

NORMALLY_HAVE = \
	unix \
	linux-64 \
	mac \
	sun

DONT_HAVE = \
	linux-32 \
	win

LIBBZ2_SYS_HAS_LIB = \
	$(addprefix $(LIBDIR)/libbz2-,$(NORMALLY_HAVE))

$(LIBBZ2_SYS_HAS_LIB):
	@ true

LIBBZ2_SYS_NEEDS_LIB = \
	$(addprefix $(LIBDIR)/libbz2-,$(DONT_HAVE))

$(LIBBZ2_SYS_NEEDS_LIB): $(LIBDIR)/libbz2


#-------------------------------------------------------------------------------
# libbz2
#
$(LIBDIR)/libbz2: $(addprefix $(LIBDIR)/libbz2.,$(LIBEXT))

LIBBZ2_SRC = \
	blocksort \
	huffman \
	crctable \
	randtable \
	compress \
	decompress \
	bzlib

LIBBZ2_OBJ = \
	$(addsuffix .$(LOBX),$(LIBBZ2_SRC))

LIBBZ2_LIB = \

$(LIBDIR)/libbz2.$(SHLX): $(LIBBZ2_OBJ)
	$(LD) --dlib --vers $(SRCDIR) -o $@ $^ $(LIBBZ2_LIB)

$(LIBDIR)/libbz2.$(LIBX): $(LIBBZ2_OBJ)
	$(LD) --slib --vers $(SRCDIR) -o $@ $^ $(LIBBZ2_LIB)


#-------------------------------------------------------------------------------
# download-$(OS)
#
download-linux: download-linux-$(BITS)

LIBBZ2_SYS_HAS_SRC = \
	$(addprefix download-,$(NORMALLY_HAVE))

$(LIBBZ2_SYS_HAS_SRC):
	@ true

LIBBZ2_SYS_NEEDS_SRC = \
	$(addprefix download-,$(DONT_HAVE))

$(LIBBZ2_SYS_NEEDS_SRC): download


#-------------------------------------------------------------------------------
# libbz2 download
#
BZIP2_VERS := \
	1.0.6

bzip2-$(BZIP2_VERS).tar.gz:
	@ echo downloading the sources
	wget http://bzip.org/$(BZIP2_VERS)/bzip2-$(BZIP2_VERS).tar.gz
	@ touch $@

libbz2.vers: bzip2-$(BZIP2_VERS).tar.gz
	@ echo unpacking the sources
	tar --no-same-owner -x -z -f $^
	mv bzip2-$(BZIP2_VERS)/*.[ch] .
	mv bzip2-$(BZIP2_VERS)/*.def .
	rm -rf bzip2-$(BZIP2_VERS)
	@ echo -ne $(BZIP2_VERS)"\n" > $@

download: libbz2.vers

.PHONY: download
