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
SCRIPT_BASE="${0%.sh}"

# os
OS="$1"
shift

# binary compiler
CC="$1"
shift

# configuration
unset TARG
unset ARGS
unset OBJX
unset SRCFILE
unset SOURCES
unset DEPENDENCIES

while [ $# -ne 0 ]
do

    case "$1" in
    --cflags)
        for ARG in $2
        do
            [ "$ARG" != "${ARG#-}" ] && ARG="/${ARG#-}"
            ARGS="$ARGS $ARG"
        done
        shift
        ;;

    --checksum)
        ;;

    --objx)
        OBJX="$2"
        shift
        ;;

    -D*)
        ARG="${1#-D}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        ARGS="$ARGS /D$ARG"
        ;;

    -I*)
        ARG="${1#-I}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        ARG="$(cygpath -w $ARG)"
        ARGS="$ARGS /I$ARG"
        ;;

    -o*)
        ARG="${1#-o}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        ARGS="$ARGS /Fo$ARG"
        TARG="${ARG%.$OBJX}"
        ;;

    -MD)
        # the /showIncludes switch will generate
        # includes to stderr, but they will need
        # to be filtered out and rewritten to the *.d
        ARGS="$ARGS /showIncludes"
        DEPENDENCIES=1
        ;;

    -fPIC|-std=c99|-ansi|-pedantic)
        ;;

    -*)
        ARG="/${1#-}"
        ARGS="$ARGS $ARG"
        ;;

    *)
        SRCFILE="$(basename $1)"
        SOURCES="$SOURCES $1"
        ARG="$(cygpath -w $1)"
        ARGS="$ARGS $ARG"
        ;;
        
    esac

    shift
done

unset STATUS

CMD="$CC $ARGS"
echo "$CMD"

# run command with redirection
if $CMD > $TARG.out 2> $TARG.err
then
    # success
    STATUS=0
else
    # failure
    STATUS=$?
fi

# check for dependencies
if [ "$DEPENDENCIES" = "1" ]
then
    sed -e '/including file/!d' -e 's/.*including file: *\([^\r\n][^\r\n]*\)/\1/g' -e '/ /d' $TARG.out > $TARG.inc
    echo -n "$TARG.$OBJX: $SOURCES" | sed -e 's/\r//g' > $TARG.d
    for inc in $(cat $TARG.inc)
    do
        echo -n " $(cygpath -u $inc)"  | sed -e 's/\r//g' >> $TARG.d
    done
fi

# repeat output files but without CR
sed -e 's/\r//g' $TARG.err > /dev/stderr
sed -e 's/\r//g' -e "/^$SRCFILE$/d" -e '/including file/d' $TARG.out

# clean up files
rm -f $TARG.out $TARG.err $TARG.inc

exit $STATUS
