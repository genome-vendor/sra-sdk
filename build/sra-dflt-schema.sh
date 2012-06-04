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

# everything except windows is fine as-is
if [ "$OS" != "win" ]
then
    echo "$CC $*"
    $CC $*
    exit $?
fi

# state
unset ARGS
unset DEPENDENCIES
unset DEPTARG
unset TARG

# process parameters for windows
while [ $# -ne 0 ]
do

    case "$1" in
    -o*)
        ARG="${1#-o}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        TARG="$ARG"
        ARG="$(cygpath -w $ARG)"
        ARGS="$ARGS -o$ARG"
        ;;

    -I*)
        ARG="${1#-I}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        ARG="$(cygpath -w $ARG)"
        ARGS="$ARGS -I$ARG"
        ;;

    -T*)
        ARG="${1#-T}"
        if [ "$ARG" = "" ]
        then
            ARG="$2"
            shift
        fi
        DEPTARG="$ARG"
        DEPENDENCIES=1
        ARG="$(cygpath -w $ARG)"
        ARGS="$ARGS -T$ARG"
        ;;

    *)
        ARG="$(cygpath -w $1)"
        ARGS="$ARGS $ARG"
        ;;

    esac

    shift

done

echo "$CC $ARGS"
if ! $CC $ARGS
then
    STATUS=$?
    rm -f $TARG $DEPTARG
    exit $?
fi

if [ $DEPENDENCIES -eq 1 ]
then
    # fix this
    rm -f $DEPTARG
fi
