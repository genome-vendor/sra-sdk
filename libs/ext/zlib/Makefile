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
MODULE = libs/ext/zlib

include $(TOP)/build/Makefile.env

INT_LIBS = \

EXT_LIBS = \
	libz

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
$(TARGDIR)/all: $(LIBDIR)/libz

.PHONY: $(TARGDIR)/all

#-------------------------------------------------------------------------------
# std
#
$(TARGDIR)/std: $(LIBDIR)/libz-$(OS)

.PHONY: $(TARGDIR)/std

#-------------------------------------------------------------------------------
# clean
#
clean: stdclean

deep-clean: clean
	rm -f *.[ch] *.def *.vers *.gz

.PHONY: clean

#-------------------------------------------------------------------------------
# libz-$(OS)
#
$(LIBDIR)/libz-linux: $(LIBDIR)/libz-linux-$(BITS)

NORMALLY_HAVE = \
	unix \
	linux-64 \
	mac \
	sun

DONT_HAVE = \
	linux-32 \
	win

LIBZ_SYS_HAS_LIB = \
	$(addprefix $(LIBDIR)/libz-,$(NORMALLY_HAVE))

$(LIBZ_SYS_HAS_LIB):
	@ true

LIBZ_SYS_NEEDS_LIB = \
	$(addprefix $(LIBDIR)/libz-,$(DONT_HAVE))

$(LIBZ_SYS_NEEDS_LIB): $(LIBDIR)/libz


#-------------------------------------------------------------------------------
# libz
#
$(LIBDIR)/libz: $(addprefix $(LIBDIR)/libz.,$(LIBEXT))

LIBZ_SRC = \
	adler32 \
	compress \
	crc32 \
	deflate \
	gzclose \
	gzlib \
	gzread \
	gzwrite \
	infback \
	inffast \
	inflate \
	inftrees \
	trees \
	uncompr \
	zutil


LIBZ_OBJ = \
	$(addsuffix .$(LOBX),$(LIBZ_SRC))

LIBZ_LIB = \

$(LIBDIR)/libz.$(SHLX): $(LIBZ_OBJ)
	$(LD) --dlib --vers $(SRCDIR) -o $@ $^ $(LIBZ_LIB)

$(LIBDIR)/libz.$(LIBX): $(LIBZ_OBJ)
	$(LD) --slib --vers $(SRCDIR) -o $@ $^ $(LIBZ_LIB)


#-------------------------------------------------------------------------------
# download-$(OS)
#
download-linux: download-linux-$(BITS)

LIBZ_SYS_HAS_SRC = \
	$(addprefix download-,$(NORMALLY_HAVE))

$(LIBZ_SYS_HAS_SRC):
	@ true

LIBZ_SYS_NEEDS_SRC = \
	$(addprefix download-,$(DONT_HAVE))

$(LIBZ_SYS_NEEDS_SRC): download


#-------------------------------------------------------------------------------
# libz download
#
ZLIB_VERS := \
	1.2.6

zlib-$(ZLIB_VERS).tar.gz:
	@ echo downloading the sources
	wget http://zlib.net/$@
	@ touch $@

libz.vers: zlib-$(ZLIB_VERS).tar.gz
	@ echo unpacking the sources
	tar --no-same-owner -x -z -f $^
	mv zlib-$(ZLIB_VERS)/*.[ch] .
	mv zlib-$(ZLIB_VERS)/win32/zlib.def libz.def
	rm -rf zlib-$(ZLIB_VERS)
	@ echo -ne $(ZLIB_VERS) "\n" > $@

download: libz.vers

.PHONY: download
