/*===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 */

#include <klib/out.h>

#include <kapp/main.h>
#include <kapp/args.h>

#include "sra-shell.vers.h"

uint32_t CC KAppVersion ( void )
{
    return SRA_SHELL_VERS;
}

rc_t CC Usage ( const Args * args )
{
    return 0;
}

const char UsageDefaultName[] = "sra-shell";

rc_t CC UsageSummary(const char *prog_name)
{
    return KOutMsg ( "Usage: %s [options] path [ path... ]\n"
                     "\n"
                     , prog_name );
}

rc_t CC KMain ( int argc, char *argv [] )
{
    return 0;
}
