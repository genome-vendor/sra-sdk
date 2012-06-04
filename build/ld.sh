#!/bin/bash
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

# prepare script name
SELF_NAME="$(basename $0)"
BUILD_DIR="$(dirname $0)"
TOP="$(dirname $BUILD_DIR)"
SCRIPT_BASE="${0%.sh}"

# os
OS="$1"
shift

# architecture
ARCH="$1"
shift

# binary loader tool
LD="$1"
shift

# configuration
unset SHLX
unset DYLX
unset LIBX
unset OBJX
unset LOBX

# parameters
TYPE=exe
STATIC=0
DYLD=0
STATICSYSLIBS=0
CHECKSUM=0
KPROC=0
THREADS=0
HAVE_GZIP=0
NEED_GZIP=0
HAVE_BZIP=0
NEED_BZIP=0
HAVE_DL=0
NEED_DL=0
unset BUILD
unset LDIRS
unset XDIRS
unset SRCDIR
unset BINDIR
unset VERSFILE
unset VERSDIR
unset TARG
unset EXT
unset OBJS
unset LIBS
unset DEPFILE

while [ $# -ne 0 ]
do

    case "$1" in
    --build)
        BUILD="$2"
        shift
        ;;

    --ldflags)
        LDFLAGS="$2"
        shift
        ;;

    --static-system-libs)
        STATICSYSLIBS=1
        ;;

    --checksum)
        CHECKSUM=1
        ;;

    --shlx)
        SHLX="$2"
        shift
        ;;

    --dylx)
        SHLX="$2"
        shift
        ;;

    --libx)
        LIBX="$2"
        shift
        ;;

    --objx)
        OBJX="$2"
        shift
        ;;

    --srcdir)
        SRCDIR="$2"
        shift
        ;;

    --bindir)
        BINDIR="$2"
        shift
        ;;

    -MD)
        DEPFILE=1
        ;;

    -L*)
        ARG="${1#-L}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        LDIRS="$LDIRS:$ARG"
        ;;
        
    -X*)
        ARG="${1#-X}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        XDIRS="$XDIRS:$ARG"
        ;;

    --dlib)
        TYPE=dlib
        ;;

    --slib)
        TYPE=slib
        ;;

    --stub)
        TYPE=stub
        ;;

    --exe)
        TYPE=exe
        ;;

    --static)
        STATIC=1
        ;;

    --vers)
        if [ -f "$2" ]
        then
            VERSFILE="$2"
        elif [ -d "$2" ]
        then
            VERSDIR="$2"
        else
            echo "$SELF_NAME: expected version file or source directory"
            exit 3
        fi
        shift
        ;;

    -o*)
        ARG="${1#-o}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        TARG="$ARG"
        ;;

    -lz|-sz|-dz)
        LIBS="$LIBS $1"
        HAVE_GZIP=1
        ;;
    -lbz2|-sbz2|-dbz2)
        LIBS="$LIBS $1"
        HAVE_BZIP=1
        ;;
    -ldl|-sdl|-ddl)
        LIBS="$LIBS $1"
        HAVE_DL=1
        ;;

    -lpthread|-spthread|-dpthread)
        THREADS=8
        ;;

    -lkfs)
        LIBS="$LIBS $1"
        if [ $STATIC -ne 0 ]
        then
            NEED_GZIP=1
            NEED_BZIP=1
            NEED_DL=1
        fi
        ;;

    -skfs)
        LIBS="$LIBS $1"
        NEED_GZIP=1
        NEED_BZIP=1
        NEED_DL=1
        ;;

    -dkfs)
        LIBS="$LIBS $1"
        NEED_GZIP=1
        NEED_BZIP=1
        NEED_DL=1
        DYLD=2
        ;;

    -lkproc)
        LIBS="$LIBS $1"
        KPROC=4
        ;;
    -skproc)
        LIBS="$LIBS $1"
        KPROC=4
        THREADS=8
        ;;
    -dkproc)
        LIBS="$LIBS $1"
        KPROC=4
        DYLD=2
        ;;

    -lncbi-bam|-sncbi-bam|-dncbi-bam)
        LIBS="$LIBS $1"
        NEED_GZIP=1
        ;;

    -l*|-s*)
        LIBS="$LIBS $1"
        ;;

    -d*)
        LIBS="$LIBS $1"
        DYLD=2
        ;;

    *.$OBJX)
        OBJS="$OBJS $1"
        ;;
        
    esac

    shift
done

# correct for prefixes
LDIRS="${LDIRS#:}"
XDIRS="${XDIRS#:}"
LIBS="${LIBS# }"
OBJS="${OBJS# }"

# split target
OUTDIR=$(dirname "$TARG")
NAME=$(basename "$TARG")

# dependency file
[ "$DEPFILE" != "" ] && DEPFILE="$NAME.$TYPE.d"

# parse target
if [ "$TYPE" = "dlib" ] && [ "$DYLX" != "" ]
then
    EXT="$DYLX"
    NAME="${NAME%.$DYLX}"
fi
unset VERS
if [ "$NAME" != "${NAME%.[0-9][0-9]*}" ]
then
    ARG="${NAME%.[0-9][0-9]*}"
    VERS="${NAME#$ARG}"
    NAME="${ARG#.}"

    if [ "$NAME" != "${NAME%.[0-9][0-9]*}" ]
    then
        ARG="${NAME%.[0-9][0-9]*}"
        VERS="${NAME#$ARG}.$VERS"
        NAME="${ARG#.}"

        if [ "$NAME" != "${NAME%.[0-9][0-9]*}" ]
        then
            ARG="${NAME%.[0-9][0-9]*}"
            VERS="${NAME#$ARG}.$VERS"
            NAME="${ARG#.}"
        fi
    fi
fi
case "$TYPE" in
dlib)
    if [ "$SHLX" != "" ]
    then
        EXT="$SHLX"
        NAME="${NAME%.$SHLX}"
    fi
    ;;
slib)
    EXT="$LIBX"
    NAME="${NAME%.$LIBX}"
esac

unset DBGAP
if [ "$NAME" != "${NAME%-dbgap}" ]
then
	DBGAP=-dbgap
	NAME="${NAME%-dbgap}"
fi

# locate version file and version
[ "$VERSDIR" != "" ] && VERSFILE="$VERSDIR/$NAME.vers"
if [ "$VERSFILE" != "" ]
then
    if [ ! -f "$VERSFILE" ]
    then
        echo "$SELF_NAME: warning - creating version file '$VERSFILE'"
        echo 1.0.0 > $VERSFILE
    fi

    if [ ! -r "$VERSFILE" ]
    then
        echo "$SELF_NAME: version file '$VERSFILE' is unreadable"
        exit 5
    fi

    ARG=$(cat $VERSFILE)
    if [ "$VERS" != "" ] && [ "$VERS" != "$ARG" ]
    then
        echo "$SELF_NAME: version from file '$VERSFILE' does not match '$VERS'"
        exit 5
    fi
    VERS="$ARG"
fi

# turn on threads for kproc
[ $STATIC -eq 0 ] && [ "$NAME" = "libkproc" ] && THREADS=8
[ $STATIC -ne 0 ] && [ $KPROC -ne 0 ] && THREADS=8

# supply missing libraries
[ $HAVE_GZIP -eq 0 ] && [ $NEED_GZIP -ne 0 ] && LIBS="$LIBS -lz"
[ $HAVE_BZIP -eq 0 ] && [ $NEED_BZIP -ne 0 ] && LIBS="$LIBS -lbz2"
[ $HAVE_DL -eq 0 ] && [ $NEED_DL -ne 0 ] && LIBS="$LIBS -ldl"

#echo "# $SELF_NAME"
#echo "#   BUILD          : $BUILD"
#echo "#   OS             : $OS"
#echo "#   ARCH           : $ARCH"
#echo "#   tool           : $LD"
#echo "#   dep file       : $DEPFILE"
#echo "#   LDFLAGS        : $LDFLAGS"
#echo "#   static sys libs: $STATICSYSLIBS"
#echo "#   checksum       : $CHECKSUM"
#echo "#   static         : $STATIC"
#echo "#   kproc          : $KPROC"
#echo "#   thread libs    : $THREADS"
#echo "#   type           : $TYPE"
#echo "#   srcdir         : $SRCDIR"
#echo "#   bindir         : $BINDIR"
#echo "#   LDIRS          : $LDIRS"
#echo "#   XDIRS          : $XDIRS"
#echo "#   vers file      : $VERSFILE"
#echo "#   vers dir       : $VERSDIR"
#echo "#   target         : $TARG"
#echo "#   outdir         : $OUTDIR"
#echo "#   name           : $NAME"
#echo "#   dbgap          : $DBGAP"
#echo "#   extension      : $EXT"
#echo "#   version        : $VERS"
#echo "#   objects        : $OBJS"
#echo "#   libraries      : $LIBS"
#echo "#   script-base    : $SCRIPT_BASE"

# overwrite dependencies
[ -f "$DEPFILE" ] && rm -f "$DEPFILE"

# generate mode
MODE=$(expr $THREADS + $KPROC + $DYLD + $STATIC)

# generate SCM flags
SCMFLAGS=$(expr $STATICSYSLIBS + $STATICSYSLIBS + $CHECKSUM)

# perform link
"$SCRIPT_BASE.$OS.$TYPE.sh" "$LD" "$ARCH" "$BUILD" "$SRCDIR" "$BINDIR" "$OUTDIR" \
    "$TARG" "$NAME" "$DBGAP" "$VERS" "$VERSFILE" "$DEPFILE" "$MODE" "$SCMFLAGS" \
    "$LDFLAGS" "$LDIRS" "$XDIRS" "$OBJS" "$LIBS" || exit $?

# establish links
if [ "$VERS" != "" ] && [ "$OS" != "win" ]
then
    $SCRIPT_BASE.$OS.ln.sh "$TYPE" "$OUTDIR" "$TARG" "$NAME" "$DBGAP" "$EXT" "$VERS"
fi

# SCM
if [ $CHECKSUM -eq 1 ] && [ "$OS" = "linux" ]
then
    # calling the scm-version-script
    # parameters are: module-name, current-md5-file, version-file
    if [ $TYPE = "dlib" ] || [ $TYPE = "exe" ] || [ $STATIC -eq 1 ]
    then
        SCM_DIR="$TOP/scm"
        LOGFILE="$SCM_DIR/scm.log"
        SCMD="$BUILD_DIR/scm.sh $NAME $TARG.md5 $VERSFILE"
        echo "$SCMD" >> $LOGFILE
        $SCMD
    fi
fi
