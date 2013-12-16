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


#include <vfs/extern.h>
#include "resolver-priv.h"

#include <vfs/manager.h>
#include <vfs/path.h>
#include <kns/manager.h>
#include <kns/http.h>
#include <kns/stream.h>
#include <kns/curl-file.h>
#include <kns/KCurlRequest.h>
#include <kfs/file.h>
#include <kfs/directory.h>
#include <kfg/repository.h>
#include <kfg/config.h>

#include <klib/text.h>
#include <klib/vector.h>
#include <klib/refcount.h>
#include <klib/namelist.h>
#include <klib/printf.h>
#include <klib/data-buffer.h>
#include <klib/debug.h> /* DBGMSG */
#include <klib/log.h>
#include <klib/rc.h>

#include <sysalloc.h>

#include <vfs/path-priv.h>
#include "path-priv.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

/* to turn off CGI name resolution for
   any refseq accessions */
#define NO_REFSEQ_CGI 1

/* to turn off CGI name resolution for
   legacy WGS packages used by refseq */
#define NO_LEGACY_WGS_REFSEQ_CGI NO_REFSEQ_CGI

#define USE_CURL 1

#define NAME_SERVICE_MAJ_VERS 1
#define NAME_SERVICE_MIN_VERS 1
#define ONE_DOT_ONE 0x01010000
#define NAME_SERVICE_VERS \
    ( ( NAME_SERVICE_MAJ_VERS << 24 ) | ( NAME_SERVICE_MIN_VERS << 16 ) )


/*--------------------------------------------------------------------------
 * String
 */
static
void CC string_whack ( void *obj, void *ignore )
{
    StringWhack ( ( String* ) obj );
}

/*--------------------------------------------------------------------------
 * VResolverAccToken
 *  breaks up an accession into parts
 *
 *  "acc" is entire accession as given
 *
 *  the remainder is divided like so:
 *
 *    [<prefix>_]<alpha><digits>[.<ext1>[.<ext2>]]
 *
 *  prefix is optional
 *  alpha can be zero length iff prefix is not zero length
 *  digits must be non-zero length
 *  ext1 and ext2 are optional
 */
typedef struct VResolverAccToken VResolverAccToken;
struct VResolverAccToken
{
    String acc;
    String prefix;
    String alpha;
    String digits;
    String ext1;
    String ext2;
    String suffix;
};

/*--------------------------------------------------------------------------
 * VResolverAlg
 *  represents a set of zero or more volumes
 *  each of which is addressed using a particular expansion algorithm
 */
typedef enum
{
    appUnknown,
    appAny,
    appREFSEQ,
    appSRA,
    appWGS,
    appNANNOT,
    appNAKMER,
    appCount
} VResolverAppID;

typedef enum
{
    algCGI,
    algLeafPath,
    algSRAFlat,
    algSRA1024,
    algSRA1000,
    algFUSE1000,
    algREFSEQ,
    algWGSFlat,
    algWGS,
    algFuseWGS,
    algSRA_NCBI,
    algSRA_EBI,

    algNANNOTFlat,
    algNANNOT,
    algFuseNANNOT,
    algNAKMERFlat,
    algNAKMER,
    algFuseNAKMER,

    /* leave as last value */
    algUnknown

} VResolverAlgID;

typedef enum
{
    cacheDisallow,
    cacheAllow
} VResolverCacheAllow;

typedef struct VResolverAlg VResolverAlg;
struct VResolverAlg
{
    /* volume paths - stored as String* */
    Vector vols;

    /* root path - borrowed reference */
    const String *root;

    /* download ticket - borrowed reference
       non-NULL means that the root is a
       resolver CGI. also, don't rely on
       presence of any volumes... */
    const String *ticket;

    /* app_id helps to filter out volumes by app */
    VResolverAppID app_id;

    /* how to expand an accession */
    VResolverAlgID alg_id;

    /* a property of the repository */
    bool protected;

    /* whether the volumes are cache-capable
       in particularl, enabled if cache forced */
    bool cache_capable;

    /* whether the volumes are cache-enabled */
    bool cache_enabled;

    /* whether the volume is disabled in config */
    bool disabled;
#if 0
    VRemoteProtocols protocols;
#endif
};


/* Whack
 */
static
void CC VResolverAlgWhack ( void *item, void *ignore )
{
    VResolverAlg *self = item;

    /* drop any volumes */
    VectorWhack ( & self -> vols, string_whack, NULL );

    /* everything else is a borrowed reference */

    free ( self );
}

/* Make
 */
static
rc_t VResolverAlgMake ( VResolverAlg **algp, const String *root,
     VResolverAppID app_id, VResolverAlgID alg_id, bool protected, bool disabled )
{
    rc_t rc;
    VResolverAlg *alg = calloc ( 1, sizeof * alg );
    if ( alg == NULL )
        rc = RC ( rcVFS, rcMgr, rcConstructing, rcMemory, rcExhausted );
    else
    {
        VectorInit ( & alg -> vols, 0, 8 );
        alg -> root = root;
        alg -> app_id = app_id;
        alg -> alg_id = alg_id;
        alg -> protected = protected;
        alg -> disabled = disabled;
        rc = 0;
    }

    * algp = alg;
    return rc;
}

/* MakeLocalWGSRefseqURI
 *  create a special URI that tells KDB how to open this
 *  obscured table, hidden away within a KAR file
 */
static
rc_t VResolverAlgMakeLocalWGSRefseqURI ( const VResolverAlg *self,
    const String *vol, const String *exp, const String *acc, const VPath ** path )
{
    if ( self -> root == NULL )
        return VPathMakeFmt ( ( VPath** ) path, NCBI_FILE_SCHEME ":%S/%S#tbl/%S", vol, exp, acc );
    return VPathMakeFmt ( ( VPath** ) path, NCBI_FILE_SCHEME ":%S/%S/%S#tbl/%S", self -> root, vol, exp, acc );
}

/* MakeeRemoteWGSRefseqURI
 *  create a special URI that tells KDB how to open this
 *  obscured table, hidden away within a KAR file
 */
static
rc_t VResolverAlgMakeRemoteWGSRefseqURI ( const VResolverAlg *self,
    const char *url, const String *acc, const VPath ** path )
{
    return VPathMakeFmt ( ( VPath** ) path, "%s#tbl/%S", url, acc );
}

/* MakeRemotePath
 *  the path is known to exist in the remote file system
 *  turn it into a VPath
 */
static
rc_t VResolverAlgMakeRemotePath ( const VResolverAlg *self,
    const char *url, const VPath ** path )
{
    return VPathMakeFmt ( ( VPath** ) path, url );
}

/* MakeLocalPath
 *  the path is known to exist in the local file system
 *  turn it into a VPath
 */
static
rc_t VResolverAlgMakeLocalPath ( const VResolverAlg *self,
    const String *vol, const String *exp, const VPath ** path )
{
    if ( self -> root == NULL )
        return VPathMakeFmt ( ( VPath** ) path, "%S/%S", vol, exp );
    return VPathMakeFmt ( ( VPath** ) path, "%S/%S/%S", self -> root, vol, exp );
}

/* expand_accession
 *  expand accession according to algorithm
 */
static
rc_t expand_algorithm ( const VResolverAlg *self, const VResolverAccToken *tok,
    char *expanded, size_t bsize, size_t *size, bool legacy_wgs_refseq )
{
    rc_t rc;
    uint32_t num;

   switch ( self -> alg_id )
    {
    case algCGI:
        return RC ( rcVFS, rcResolver, rcResolving, rcType, rcIncorrect );
    case algLeafPath:
        rc = string_printf ( expanded, bsize, size, "%S", & tok -> acc );
        break;
    case algSRAFlat:
        rc = string_printf ( expanded, bsize, size,
            "%S%S.sra", & tok -> alpha, & tok -> digits );
        break;
    case algSRA1024:
        num = ( uint32_t ) strtoul ( tok -> digits . addr, NULL, 10 );
        rc = string_printf ( expanded, bsize, size,
            "%S/%06u/%S%S.sra", & tok -> alpha, num >> 10, & tok -> alpha, & tok -> digits );
        break;
    case algSRA1000:
        num = ( uint32_t ) ( tok -> alpha . size + tok -> digits . size - 3 );
        rc = string_printf ( expanded, bsize, size,
            "%S/%.*S/%S%S.sra", & tok -> alpha, num, & tok -> acc, & tok -> alpha, & tok -> digits );
        break;
    case algFUSE1000:
        num = ( uint32_t ) ( tok -> alpha . size + tok -> digits . size - 3 );
        rc = string_printf ( expanded, bsize, size,
            "%S/%.*S/%S%S/%S%S.sra", & tok -> alpha, num, & tok -> acc, 
            & tok -> alpha, & tok -> digits, & tok -> alpha, & tok -> digits );
        break;
    case algREFSEQ:
        if ( ! legacy_wgs_refseq )
            rc = string_printf ( expanded, bsize, size, "%S", & tok -> acc );
        else
            rc = string_printf ( expanded, bsize, size, "%S%.2S", & tok -> alpha, & tok -> digits );
        break;
    case algWGSFlat:
        num = ( uint32_t ) ( tok -> alpha . size + 2 );
        if ( tok -> prefix . size != 0 )
            num += tok -> prefix . size + 1;
        rc = string_printf ( expanded, bsize, size,
            "%.*S", num, & tok -> acc );
        break;
    case algWGS:
        num = ( uint32_t ) ( tok -> alpha . size + 2 );
        if ( tok -> prefix . size != 0 )
            num += tok -> prefix . size + 1;
        rc = string_printf ( expanded, bsize, size,
            "WGS/%.2s/%.2s/%.*S", tok -> alpha . addr, tok -> alpha . addr + 2, num, & tok -> acc );
        break;
    case algFuseWGS:
        num = ( uint32_t ) ( tok -> alpha . size + 2 );
        if ( tok -> prefix . size != 0 )
            num += tok -> prefix . size + 1;
        rc = string_printf ( expanded, bsize, size,
            "%.2s/%.2s/%.*S", tok -> alpha . addr, tok -> alpha . addr + 2, num, & tok -> acc );
        break;
    case algSRA_NCBI:
        num = ( uint32_t ) strtoul ( tok -> digits . addr, NULL, 10 );
        rc = string_printf ( expanded, bsize, size,
            "%S/%06u/%S%S", & tok -> alpha, num >> 10, & tok -> alpha, & tok -> digits );
        break;
    case algSRA_EBI:
        num = ( uint32_t ) ( tok -> alpha . size + tok -> digits . size - 3 );
        rc = string_printf ( expanded, bsize, size,
            "%S/%.*S/%S%S", & tok -> alpha, num, & tok -> acc, & tok -> alpha, & tok -> digits );
        break;

    case algNANNOTFlat:
        rc = string_printf ( expanded, bsize, size, "%S", & tok -> acc );
        break;
    case algNANNOT:
        num = ( uint32_t ) strtoul ( tok -> digits . addr, NULL, 10 );
        rc = string_printf ( expanded, bsize, size,
            "%03u/%03u/%S", num / 1000000, ( num / 1000 ) % 1000, & tok -> acc );
        break;
    case algFuseNANNOT:
        num = ( uint32_t ) strtoul ( tok -> digits . addr, NULL, 10 );
        rc = string_printf ( expanded, bsize, size,
            "%03u/%03u/%S", num / 1000000, ( num / 1000 ) % 1000, & tok -> acc );
        break;

    case algNAKMERFlat:
        rc = string_printf ( expanded, bsize, size, "%S", & tok -> acc );
        break;
    case algNAKMER:
        num = ( uint32_t ) strtoul ( tok -> digits . addr, NULL, 10 );
        rc = string_printf ( expanded, bsize, size,
            "kmer/%03u/%03u/%S", num / 1000000, ( num / 1000 ) % 1000, & tok -> acc );
        break;
    case algFuseNAKMER:
        num = ( uint32_t ) strtoul ( tok -> digits . addr, NULL, 10 );
        rc = string_printf ( expanded, bsize, size,
            "kmer/%03u/%03u/%S", num / 1000000, ( num / 1000 ) % 1000, & tok -> acc );
        break;

    default:
        return RC ( rcVFS, rcResolver, rcResolving, rcType, rcUnrecognized );
    }

   return rc;
}

/* LocalResolve
 *  resolve an accession into a VPath or not found
 *
 *  1. expand accession according to algorithm
 *  2. search all volumes for accession
 *  3. return not found or new VPath
 */
static
rc_t VResolverAlgLocalResolve ( const VResolverAlg *self,
    const KDirectory *wd, const VResolverAccToken *tok,
    const VPath ** path, bool legacy_wgs_refseq, bool for_cache )
{
    KPathType kpt;
    uint32_t i, count;

    /* expanded accession */
    String exp;
    size_t size;
    char expanded [ 256 ];

    /* in some cases, "root" is NULL */
    const String *vol, *root = self -> root;

    /* expand the accession */
    rc_t rc = expand_algorithm ( self, tok, expanded, sizeof expanded, & size, legacy_wgs_refseq );

    /* should never have a problem here... */
    if ( rc != 0 )
        return rc;

    /* if this is to detect a cache file, append extension */
    if ( for_cache )
    {
        size += string_copy ( & expanded [ size ], sizeof expanded - size, ".cache", sizeof ".cache" - 1 );
        if ( size == sizeof expanded )
            return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    }

    /* turn the expanded portion into a String
       we know that size is also length due to
       accession content rules */
    StringInit ( & exp, expanded, size, ( uint32_t ) size );

    /* remove the cache extension */
    if ( for_cache )
    {
        exp . len -= sizeof ".cache" - 1;
        exp . size -= sizeof ".cache" - 1;
    }

    /* now search all volumes */
    count = VectorLength ( & self -> vols );
    if ( root == NULL )
    {
        for ( i = 0; i < count; ++ i )
        {
            vol = VectorGet ( & self -> vols, i );
            kpt = KDirectoryPathType ( wd, "%.*s/%.*s"
                , ( int ) vol -> size, vol -> addr
                , ( int ) size, expanded );
            switch ( kpt & ~ kptAlias )
            {
            case kptFile:
            case kptDir:
                if ( legacy_wgs_refseq )
                    return VResolverAlgMakeLocalWGSRefseqURI ( self, vol, & exp, & tok -> acc, path );
                return VResolverAlgMakeLocalPath ( self, vol, & exp, path );
            default:
                break;
            }
        }
    }
    else
    {
        for ( i = 0; i < count; ++ i )
        {
            vol = VectorGet ( & self -> vols, i );
            kpt = KDirectoryPathType ( wd, "%.*s/%.*s/%.*s"
                , ( int ) root -> size, root -> addr
                , ( int ) vol -> size, vol -> addr
                , ( int ) size, expanded );
            switch ( kpt & ~ kptAlias )
            {
            case kptFile:
            case kptDir:
                if ( legacy_wgs_refseq )
                    return VResolverAlgMakeLocalWGSRefseqURI ( self, vol, & exp, & tok -> acc, path );
                return VResolverAlgMakeLocalPath ( self, vol, & exp, path );
            default:
                break;
            }
        }
    }
    
    return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}

static
rc_t VPathCheckFromNamesCGI ( const VPath * path, const String *ticket, const VPath ** mapping )
{
    size_t i, size;
    const char * start;

    /* must have an explicit scheme */
    if ( ! path -> from_uri )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

    /* can only be http or fasp */
    switch ( path -> scheme_type )
    {
    case vpuri_http:
    case vpuri_fasp:
        break;
    default:
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    }


    /* must have a host-spec with all ascii-characters */
    switch ( path -> host_type )
    {
    case vhDNSName:
        if ( path -> host . size == 0 || path -> host . size != ( size_t ) path -> host . len )
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
        start = path -> host . addr;
        size = path -> host . size;
        for ( i = 0; i < size; ++ i )
        {
            if ( isalnum ( start [ i ] ) )
                continue;
            switch ( start [ i ] )
            {
            case '.':
            case '-':
                continue;
            }
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
        }
        break;
    case vhIPv4:
    case vhIPv6:
        break;
    }

    /* must have a full-path */
    if ( path -> path_type != vpFullPath )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    /* only ascii characters */
    assert ( path -> path . size != 0 );
    if ( path -> path . size != ( size_t ) path -> path . len )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    start = path -> path . addr;
    size = path -> path . size;
    for ( i = 0; i < size; ++ i )
    {
        if ( isalnum ( start [ i ] ) )
            continue;
        switch ( start [ i ] )
        {
        case '/':
        case '.':
        case '-':
            continue;
        }
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    }

#if DO_NOT_USE_TIC_HACK
    /* if the ticket was placed into the mapped path */
    if ( mapping != NULL )
        ticket = NULL;
#endif

    if ( path -> query . size != 0 )
    {
        String name, val, req;

        /* query must match ticket */
        if ( ticket == NULL )
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

        StringSubstr ( & path -> query, & name, 0, 5 );
        StringSubstr ( & path -> query, & val, 5, 0 );
        if ( ! StringEqual ( & val, ticket ) )
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
        CONST_STRING ( & req, "?tic=" );
        if ( ! StringEqual ( & name, & req ) )
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    }

    /* cannot have a fragment */
    if ( path -> fragment . size != 0 )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

    return 0;
}


/* ParseResolverCGIResponse_1_0
 *  expect single row table, with this structure:
 *
 *  <accession>|<download-ticket>|<url>|<result-code>|<message>
 */
static
rc_t VResolverAlgParseResolverCGIResponse_1_0 ( const char *start, size_t size,
    const VPath ** path, const String *acc, const String *ticket )
{
    rc_t rc;
    KLogLevel lvl;
    char *rslt_end;
    uint32_t result_code;

    String accession, download_ticket, url, rslt_code, msg;

    /* get accession */
    const char *end = start + size;
    const char *sep = string_chr ( start, size, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & accession, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get download-ticket */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & download_ticket, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get url */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & url, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get result-code */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & rslt_code, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get msg */
    start = sep + 1;
    for ( sep = end; sep > start; -- sep )
    {
        switch ( sep [ -1 ] )
        {
        case '\n':
        case '\r':
            continue;
        default:
            break;
        }

        break;
    }
    StringInit ( & msg, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* compare acc to accession */
    if ( ! StringEqual ( & accession, acc ) )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

    /* compare ticket
       currently this makes sense with 1 request from a known workspace */
    if ( download_ticket . size != 0 )
    {
        if ( ticket == NULL || ! StringEqual ( & download_ticket, ticket ) )
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    }

    /* get the result code */
    if ( rslt_code . size == 0 )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    result_code = strtoul ( rslt_code . addr, & rslt_end, 10 );
    if ( ( const char* ) rslt_end - rslt_code . addr != rslt_code . size )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

    /* still have to test the URL */    

    switch ( result_code / 100 )
    {
    case 1:
        /* informational response
           not much we can do here */
        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        break;

    case 2:
        /* successful response
           but can only handle 200 */
        if ( result_code == 200 )
        {
            /* normal public response */
            if ( download_ticket . size == 0 )
                rc = VPathMakeFmt ( ( VPath** ) path, "%S", & url );
            else
            {
                /* protected response */
                rc = VPathMakeFmt ( ( VPath** ) path, "%S?tic=%S", & url, & download_ticket );
            }

            if ( rc == 0 )
            {
                rc = VPathCheckFromNamesCGI ( * path, ticket, NULL );
                if ( rc == 0 )
                    return 0;

                VPathRelease ( * path );
                * path = NULL;
            }

            return rc;
        }

        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        break;

    case 3:
        /* redirection
           currently this is being handled by our request object */
        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        break;

    case 4:
        /* client error */
        lvl = klogErr;
        switch ( result_code )
        {
        case 400:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcInvalid );
            break;
        case 401:
        case 403:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcQuery, rcUnauthorized );
            break;
        case 404:
            return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
        case 410:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
            break;
        default:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        }
        break;

    case 5:
        /* server error */
        lvl = klogSys;
        switch ( result_code )
        {
        case 503:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcDatabase, rcNotAvailable );
            break;
        case 504:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcTimeout, rcExhausted );
            break;
        default:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        }
        break;

    default:
        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
    }

    /* log message to user */
    PLOGERR ( lvl, ( lvl, rc, "failed to resolve accession '$(acc)' - $(msg) ( $(code) )",
        "acc=%S,msg=%S,code=%u", acc, & msg, result_code ) );
    return rc;
}


/* ParseResolverCGIResponse_1_1
 *  expect single row table, with this structure (SRA-1690) :
 *
 *  <accession>|obj-id|name|size|mod-date|md5|<download-ticket>|<url>|<result-code>|<message>
 */
static
rc_t VResolverAlgParseResolverCGIResponse_1_1 ( const char *start, size_t size,
    const VPath ** path, const VPath ** mapping, const String *acc, const String *ticket )
{
    rc_t rc;
    KLogLevel lvl;
    char *rslt_end;
    uint32_t result_code;

    String accession, obj_id, name, size_str, mod_date, md5, download_ticket, url, rslt_code, msg;

    /* get accession */
    const char *end = start + size;
    const char *sep = string_chr ( start, size, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & accession, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get obj-id */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & obj_id, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get name */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & name, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get size */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & size_str, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get mod-date */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & mod_date, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get md5 */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & md5, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get download-ticket */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & download_ticket, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get url */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & url, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get result-code */
    start = sep + 1;
    sep = string_chr ( start, end - start, '|' );
    if ( sep == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
    StringInit ( & rslt_code, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* get msg */
    start = sep + 1;
    for ( sep = end; sep > start; -- sep )
    {
        switch ( sep [ -1 ] )
        {
        case '\n':
        case '\r':
            continue;
        default:
            break;
        }

        break;
    }
    StringInit ( & msg, start, sep - start, ( uint32_t ) ( sep - start ) );

    /* compare acc to accession or obj_id */
    if ( ! StringEqual ( & accession, acc ) && ! StringEqual ( & obj_id, acc ) )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

    /* compare ticket
       currently this makes sense with 1 request from a known workspace */
    if ( download_ticket . size != 0 )
    {
        if ( ticket == NULL || ! StringEqual ( & download_ticket, ticket ) )
            return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    }

    /* get the result code */
    if ( rslt_code . size == 0 )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );
    result_code = strtoul ( rslt_code . addr, & rslt_end, 10 );
    if ( ( const char* ) rslt_end - rslt_code . addr != rslt_code . size )
        return RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcCorrupt );

    /* still have to test the URL */    

    switch ( result_code / 100 )
    {
    case 1:
        /* informational response
           not much we can do here */
        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        break;

    case 2:
        /* successful response
           but can only handle 200 */
        if ( result_code == 200 )
        {
            /* normal public response */
            if ( download_ticket . size == 0
#if DO_NOT_USE_TIC_HACK
                 || mapping != NULL
#endif
                )
            {
                rc = VPathMakeFmt ( ( VPath** ) path, "%S", & url );
            }
            else
            {
                /* protected response */
                rc = VPathMakeFmt ( ( VPath** ) path, "%S?tic=%S", & url, & download_ticket );
            }

            if ( rc == 0 )
            {
                rc = VPathCheckFromNamesCGI ( * path, ticket, mapping );
                if ( rc == 0 )
                {
                    if ( mapping == NULL )
                        return 0;

                    if ( download_ticket . size != 0 )
                    {
                        if ( accession . size != 0 )
                            rc = VPathMakeFmt ( ( VPath** ) mapping, "ncbi-acc:%S?tic=%S", & accession, & download_ticket );
                        else if ( name . size == 0 )
                            return 0;
                        else
                            rc = VPathMakeFmt ( ( VPath** ) mapping, "ncbi-file:%S?tic=%S", & name, & download_ticket );
                    }
                    else if ( accession . size != 0 )
                        rc = VPathMakeFmt ( ( VPath** ) mapping, "ncbi-acc:%S", & accession );
                    else if ( name . size == 0 )
                        return 0;
                    else
                        rc = VPathMakeFmt ( ( VPath** ) mapping, "ncbi-file:%S", & name );

                    if ( rc == 0 )
                        return 0;
                }

                VPathRelease ( * path );
                * path = NULL;
            }

            return rc;
        }

        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        break;

    case 3:
        /* redirection
           currently this is being handled by our request object */
        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        break;

    case 4:
        /* client error */
        lvl = klogErr;
        switch ( result_code )
        {
        case 400:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcMessage, rcInvalid );
            break;
        case 401:
        case 403:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcQuery, rcUnauthorized );
            break;
        case 404:
            return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
        case 410:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
            break;
        default:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        }
        break;

    case 5:
        /* server error */
        lvl = klogSys;
        switch ( result_code )
        {
        case 503:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcDatabase, rcNotAvailable );
            break;
        case 504:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcTimeout, rcExhausted );
            break;
        default:
            rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
        }
        break;

    default:
        lvl = klogInt;
        rc = RC ( rcVFS, rcResolver, rcResolving, rcError, rcUnexpected );
    }

    /* log message to user */
    PLOGERR ( lvl, ( lvl, rc, "failed to resolve accession '$(acc)' - $(msg) ( $(code) )",
        "acc=%S,msg=%S,code=%u", acc, & msg, result_code ) );
    return rc;
}


/* ParseResolverCGIResponse
 *  the response should be NUL terminated
 *  but should also be close to the size of result
 */
static
rc_t VResolverAlgParseResolverCGIResponse ( const KDataBuffer *result,
    const VPath ** path, const VPath ** mapping, const String *acc, const String *ticket )
{
    /* the textual response */
    const char *start = ( const void* ) result -> base;
    size_t i, size = KDataBufferBytes ( result );

    DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS), (" Response = %s\n", start));

    /* peel back buffer to significant bytes */
    while ( size > 0 && start [ size - 1 ] == 0 )
        -- size;

    /* skip over blanks */
    for ( i = 0; i < size; ++ i )
    {
        if ( ! isspace ( start [ i ] ) )
            break;
    }

    /* at this point, we expect only version 1.0 ... */
    if ( string_cmp ( & start [ i ], size - i, "#1.0", sizeof "#1.0" - 1, sizeof "#1.0" - 1 ) == 0 )
    {
        do
        {
            /* accept version line */
            i += sizeof "#1.0" - 1;

            /* must be followed by eoln */
            if ( start [ i ] == '\r' && start [ i + 1 ] == '\n' )
                i += 2;
            else if ( start [ i ] == '\n' )
                i += 1;
            else
                break;

            /* parse 1.0 response table */
            return VResolverAlgParseResolverCGIResponse_1_0 ( & start [ i ], size - i, path, acc, ticket );
        }
        while ( false );
    }

    /* ... and 1.1 */
    if ( string_cmp ( & start [ i ], size - i, "#1.1", sizeof "#1.1" - 1, sizeof "#1.1" - 1 ) == 0 )
    {
        do
        {
            /* accept version line */
            i += sizeof "#1.1" - 1;

            /* must be followed by eoln */
            if ( start [ i ] == '\r' && start [ i + 1 ] == '\n' )
                i += 2;
            else if ( start [ i ] == '\n' )
                i += 1;
            else
                break;

            /* parse 1.0 response table */
            return VResolverAlgParseResolverCGIResponse_1_1 ( & start [ i ], size - i, path, mapping, acc, ticket );
        }
        while ( false );
    }

    return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}

/* RemoteProtectedResolve
 *  use NCBI CGI to resolve accession into URL
 */
static
rc_t VResolverAlgRemoteProtectedResolve ( const VResolverAlg *self,
    const KNSManager *kns, VRemoteProtocols protocols, const String *acc,
    const VPath ** path, const VPath ** mapping, bool legacy_wgs_refseq )
{
#if USE_CURL
    struct KCurlRequest *req;
    rc_t rc = 0;
    DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS), ("names.cgi = %S\n", self -> root));
    rc = KNSManagerMakeCurlRequest ( kns, & req, self -> root -> addr, false );
    if ( rc == 0 )
    {
        String name, val;

        /* build up POST information: */
        CONST_STRING ( & name, "version" );
#if NAME_SERVICE_VERS == ONE_DOT_ONE
        CONST_STRING (& val, "1.1" );
#else
        CONST_STRING (& val, "1.0" );
#endif
        DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS), ("  %S = %S\n", &name, &val));
        rc = KCurlRequestAddSField ( req, & name, & val );
        if ( rc == 0 )
        {
            CONST_STRING ( & name, "acc" );
            DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS), ("  %S = %S\n", &name, acc));
            rc = KCurlRequestAddSField ( req, & name, acc );
        }
        if ( rc == 0 && legacy_wgs_refseq )
        {
            CONST_STRING ( & name, "ctx" );
            CONST_STRING (& val, "refseq" );
            DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS), ("  %S = %S\n", &name, &val));
            rc = KCurlRequestAddSField ( req, & name, & val );
        }
        if ( rc == 0 && self -> ticket != NULL )
        {
            CONST_STRING ( & name, "tic" );
            DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS),
                ("  %S = %S\n", &name, self->ticket));
            rc = KCurlRequestAddSField ( req, & name, self -> ticket );
        }
#if NAME_SERVICE_VERS >= ONE_DOT_ONE /* SRA-1690 */
        if ( rc == 0 )
        {
            CONST_STRING ( & name, "accept-proto" );
            switch ( protocols )
            {
            case eProtocolHttp:
                CONST_STRING ( & val, "http" );
                break;
            case eProtocolFasp:
                CONST_STRING ( & val, "fasp" );
                break;
            case eProtocolFaspHttp:
                CONST_STRING ( & val, "fasp,http" );
                break;
            case eProtocolHttpFasp:
                CONST_STRING ( & val, "http,fasp" );
                break;
            default:
                rc = RC ( rcVFS, rcResolver, rcResolving, rcParam, rcInvalid );
            }
            if ( rc == 0 ) {
                DBGMSG(DBG_VFS, DBG_FLAG(DBG_VFS),
                    ("  %S = %S\n", &name, &val));
                rc = KCurlRequestAddSField ( req, & name, & val );
            }
        }
#endif

        /* execute post */
        if ( rc == 0 )
        {
            KDataBuffer result;
            memset ( & result, 0, sizeof result );
            rc = KCurlRequestPerform ( req, & result );
            if ( rc == 0 )
            {
                /* expect a table as a NUL-terminated string, but
                   having close to the number of bytes in results */
                rc = VResolverAlgParseResolverCGIResponse ( & result, path, mapping, acc, self -> ticket );
                KDataBufferWhack ( & result );
            }
        }

        KCurlRequestRelease ( req );
    }
#else
    KHttpRequest *req;
    rc_t rc = KNSManagerMakeRequest ( kns, &req, 0x01000000, NULL, self -> root -> addr );
    if ( rc == 0 )
    {
        /* build up POST information: */
        rc = KHttpRequestAddPostParam ( req, "version=%u.%u",
            NAME_SERVICE_MAJ_VERS, NAME_SERVICE_MIN_VERS );
        if ( rc == 0 )
            rc = KHttpRequestAddPostParam ( req, "acc=%S", acc ); 
        if ( rc == 0 && legacy_wgs_refseq )
            rc = KHttpRequestAddPostParam ( req, "ctx=refseq" );
        if ( rc == 0 && self -> ticket != NULL )
            rc = KHttpRequestAddPostParam ( req, "tic=%S", self -> ticket );
        if ( NAME_SERVICE_VERS >= ONE_DOT_ONE )
        {
            const char *val;
            switch ( protocols )
            {
            case eProtocolHttp:
                val = "http";
                break;
            case eProtocolFasp:
                val = "fasp";
                break;
            case eProtocolFaspHttp:
                val = "fasp,http";
                break;
            case eProtocolHttpFasp:
                val = "http,fasp";
                break;
            default:
                val = NULL;
                rc = RC ( rcVFS, rcResolver, rcResolving, rcParam, rcInvalid );
            }

            if ( rc == 0 )
                rc = KHttpRequestAddPostParam ( req, "accept-proto=%s", val );
        }

        if ( rc == 0 )
        {
            KHttpResult *rslt;
            
            rc = KHttpRequestPOST ( req, &rslt );
            if ( rc == 0 )
            {
                uint32_t code;

                rc = KHttpResultStatus ( rslt, &code, NULL, 0, NULL );
                if ( code == 200 )
                {
                    KStream *response;
                    
                    rc = KHttpResultGetInputStream ( rslt, &response );
                    if ( rc == 0 )
                    {
                        KDataBuffer result;
                        size_t num_read;
                        size_t total = 0;
                        
                        memset ( &result, 0, sizeof result );
                        KDataBufferMakeBytes ( & result, 4096 );

                        while ( 1 )
                        {
                            uint8_t *base;
                            uint64_t avail = result . elem_count - total;
                            if ( avail < 256 )
                            {
                                rc = KDataBufferResize ( & result, result . elem_count + 4096 );
                                if ( rc != 0 )
                                    break;
                            }
                            
                            base = result . base;
                            rc = KStreamRead ( response, & base [ total ], result . elem_count - total, &num_read );
                            if ( rc != 0 )
                            {
                                /* TBD - look more closely at rc */
                                if ( num_read > 0 )
                                    rc = 0;
                                else
                                    break;
                            }
                            if ( num_read == 0 )
                                break;
                        }
                        
                        if ( rc == 0 )
                        {
                            rc = VResolverAlgParseResolverCGIResponse ( & result, path, mapping, acc, self -> ticket );
                            KDataBufferWhack ( &result );
                        }

                        KStreamRelease ( response );
                    }
                }
                KHttpResultRelease ( rslt );
            }
        }
        KHttpRequestRelease ( req );
    }
#endif
    return rc == 0 ? 0 : RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}

/* RemoteResolve
 *  resolve an accession into a VPath or not found
 *
 *  1. expand accession according to algorithm
 *  2. search all volumes for accession
 *  3. return not found or new VPath
 */
static
rc_t VResolverAlgRemoteResolve ( const VResolverAlg *self,
    const KNSManager *kns, VRemoteProtocols protocols, const VResolverAccToken *tok,
    const VPath ** path, const VPath ** mapping, const KFile ** opt_file_rtn, bool legacy_wgs_refseq )
{
    rc_t rc;
    uint32_t i, count;

    /* expanded accession */
    String exp;
    size_t size;
    char expanded [ 256 ];

    const String *root;

    /* check for download ticket */
    if ( self -> alg_id == algCGI
#if NO_LEGACY_WGS_REFSEQ_CGI
         && ! legacy_wgs_refseq
#endif
        )
    {
        rc = VResolverAlgRemoteProtectedResolve ( self,
            kns, protocols, & tok -> acc, path, mapping, legacy_wgs_refseq );
        if (rc == 0 && path != NULL && *path != NULL &&
            opt_file_rtn != NULL && *opt_file_rtn == NULL &&
            VPathGetUri_t ( * path ) != vpuri_fasp )
        {
            const String *s = NULL;
            rc_t rc = VPathMakeString(*path, &s);
            if (rc != 0) {
                LOGERR(klogInt, rc,
                    "failed to make string from remote protected path");
            }
            else {
#if USE_CURL
                rc = KCurlFileMake(opt_file_rtn, s->addr, false);
#else
                rc = KNSManagerMakeHttpFile ( kns, opt_file_rtn, NULL, 0x01010000, "%S", s );
#endif
                if (rc != 0) {
                    PLOGERR(klogInt, (klogInt, rc,
                        "failed to open file for $(path)", "path=%s", s->addr));
                }
                free((void*)s);
            }
        }
        return rc;
    }

    /* for remote, root can never be NULL */
    root = self -> root;

    /* expand the accession */
    rc = expand_algorithm ( self, tok, expanded, sizeof expanded, & size, legacy_wgs_refseq );

    /* should never have a problem here... */
    if ( rc != 0 )
        return rc;

    /* turn the expanded portion into a String
       we know that size is also length due to
       accession content rules */
    StringInit ( & exp, expanded, size, ( uint32_t ) size );

    /* now search all remote volumes */
    count = VectorLength ( & self -> vols );
    for ( i = 0; i < count; ++ i )
    {
        char url [ 8192 ];
        const String *vol = VectorGet ( & self -> vols, i );
        rc = string_printf ( url, sizeof url, NULL, "%S/%S/%S", root, vol, & exp );
        if ( rc == 0 )
        {
            const KFile *f;
#if USE_CURL
            rc = KCurlFileMake ( & f, url, false );
#else
            rc = KNSManagerMakeHttpFile ( kns, & f, NULL, 0x01010000, url );
#endif
            if ( rc == 0 )
            {
                if ( opt_file_rtn != NULL )
                    * opt_file_rtn = f;
                else
                    KFileRelease ( f );

                if ( legacy_wgs_refseq )
                    return VResolverAlgMakeRemoteWGSRefseqURI ( self, url, & tok -> acc, path );
                return VResolverAlgMakeRemotePath ( self, url, path );
            }
        }
    }
    
    return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}


/* CacheResolve
 *  try to resolve accession for currently cached file
 */
static
rc_t VResolverAlgCacheResolve ( const VResolverAlg *self,
    const KDirectory *wd, const VResolverAccToken *tok,
    const VPath ** path, bool legacy_wgs_refseq )
{
    /* see if the cache file already exists */
    const bool for_cache = true;
    rc_t rc = VResolverAlgLocalResolve ( self, wd, tok, path, legacy_wgs_refseq, for_cache );
    if ( rc == 0 )
        return 0;

    /* TBD - see if any of these volumes is a good candidate for creating a file */

    return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}



/* MakeCachePath
 *  we have an accession that matches this volume
 *  create a path for it
 */
static
rc_t VResolverAlgMakeCachePath ( const VResolverAlg *self,
    const VResolverAccToken *tok, const VPath ** path, bool legacy_wgs_refseq )
{
    uint32_t i, count;

    /* expanded accession */
    String exp;
    size_t size;
    char expanded [ 256 ];

    const String *vol;

    /* expand the accession */
    rc_t rc = expand_algorithm ( self, tok, expanded, sizeof expanded, & size, legacy_wgs_refseq );

    /* should never have a problem here... */
    if ( rc != 0 )
        return rc;

    /* turn the expanded portion into a String
       we know that size is also length due to
       accession content rules */
    StringInit ( & exp, expanded, size, ( uint32_t ) size );

    /* now search all volumes */
    count = VectorLength ( & self -> vols );
    for ( i = 0; i < count; ++ i )
    {
        vol = VectorGet ( & self -> vols, i );
        return VResolverAlgMakeLocalPath ( self, vol, & exp, path );
    }
    
    return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );
}



/*--------------------------------------------------------------------------
 * VResolver
 */
struct VResolver
{
    /* root paths - stored as String* */
    Vector roots;

    /* volume algorithms - stored as VResolverAlg* */
    Vector local;
    Vector remote;

    /* working directory for testing local paths */
    const KDirectory *wd;

    /* if there is a protected remote resolver,
       we will need a KNS manager */
    const KNSManager *kns;

    /* if there is a working protected repository,
       store the download ticket here */
    const String *ticket;

    KRefcount refcount;

    /* counters for various app volumes */
    uint32_t num_app_vols [ appCount ];

    /* preferred protocols preferences. Default: HTTP */
    VRemoteProtocols protocols;
};


/* "process" global settings
 *  actually, they are library-closure global
 */
static atomic32_t enable_local, enable_remote, enable_cache;


/* Whack
 */
static
rc_t VResolverWhack ( VResolver *self )
{
    KRefcountWhack ( & self -> refcount, "VResolver" );

    /* drop all remote volume sets */
    VectorWhack ( & self -> remote, VResolverAlgWhack, NULL );

    /* drop local volume sets */
    VectorWhack ( & self -> local, VResolverAlgWhack, NULL );

    /* drop download ticket */
    if ( self -> ticket != NULL )
        StringWhack ( ( String* ) self -> ticket );

    /* drop root paths */
    VectorWhack ( & self -> roots, string_whack, NULL );

    /* release kns */
    if ( self -> kns != NULL )
        KNSManagerRelease ( self -> kns );

    /* release directory onto local file system */
    KDirectoryRelease ( self -> wd );

    free ( self );
    return 0;
}


/* AddRef
 * Release
 */
LIB_EXPORT
rc_t CC VResolverAddRef ( const VResolver * self )
{
    if ( self != NULL )
    {
        switch ( KRefcountAdd ( & self -> refcount, "VResolver" ) )
        {
        case krefOkay:
            break;
        case krefZero:
            return RC ( rcVFS, rcResolver, rcAttaching, rcRefcount, rcIncorrect );
        case krefLimit:
            return RC ( rcVFS, rcResolver, rcAttaching, rcRefcount, rcExhausted );
        case krefNegative:
            return RC ( rcVFS, rcResolver, rcAttaching, rcRefcount, rcInvalid );
        default:
            return RC ( rcVFS, rcResolver, rcAttaching, rcRefcount, rcUnknown );
        }
    }
    return 0;
}

LIB_EXPORT
rc_t CC VResolverRelease ( const VResolver * self )
{
    rc_t rc = 0;
    if ( self != NULL )
    {
        switch ( KRefcountDrop ( & self -> refcount, "VResolver" ) )
        {
        case krefOkay:
        case krefZero:
            break;
        case krefWhack:
            VResolverWhack ( ( VResolver* ) self );
            break;
        case krefNegative:
            return RC ( rcVFS, rcResolver, rcAttaching, rcRefcount, rcInvalid );
        default:
            rc = RC ( rcVFS, rcResolver, rcAttaching, rcRefcount, rcUnknown );
            break;            
        }
    }
    return rc;
}


/* get_accession_code
 */
static
uint32_t get_accession_code ( const String * accession, VResolverAccToken *tok )
{
#if USE_VPATH_OPTIONS_STRINGS
#error this thing is wrong
#else

#define MAX_ACCESSION_LEN 20

    uint32_t code;

    const char *acc = accession -> addr;
    size_t i, size = accession -> size;

    /* capture the whole accession */
    tok -> acc = * accession;

    /* scan prefix or alpha */
    for ( i = 0; i < size; ++ i )
    {
        if ( ! isalpha ( acc [ i ] ) )
            break;
    }

    /* terrible situation - unrecognizable */
    if ( i == size || i == 0 || i >= MAX_ACCESSION_LEN )
    {
        StringInit ( & tok -> prefix, acc, 0, 0 );
        StringInit ( & tok -> alpha, acc, i, i );
        StringInit ( & tok -> digits, & acc [ i ], 0, 0 );
        tok -> ext1 = tok -> ext2 = tok -> suffix = tok -> digits;
        return 0;
    }

    /* if stopped on '_', we have a prefix */
    if ( acc [ i ] == '_' )
    {
        /* prefix
           store only its presence, not length */
        code = 1 << 4 * 4;
        StringInit ( & tok -> prefix, acc, i, i );

        /* remove prefix */
        acc += ++ i;
        size -= i;

        /* scan for alpha */
        for ( i = 0; i < size; ++ i )
        {
            if ( ! isalpha ( acc [ i ] ) )
                break;
        }

        if ( i == size || i >= MAX_ACCESSION_LEN )
        {
            StringInit ( & tok -> alpha, acc, i, i );
            StringInit ( & tok -> digits, & acc [ i ], 0, 0 );
            tok -> ext1 = tok -> ext2 = tok -> suffix = tok -> digits;
            return 0;
        }

        code |= i << 4 * 3;
        StringInit ( & tok -> alpha, acc, i, i );
    }
    else if ( ! isdigit ( acc [ i ] ) )
    {
        StringInit ( & tok -> prefix, acc, 0, 0 );
        StringInit ( & tok -> alpha, acc, i, i );
        StringInit ( & tok -> digits, & acc [ i ], 0, 0 );
        tok -> ext1 = tok -> ext2 = tok -> suffix = tok -> digits;
        return 0;
    }
    else
    {
        /* alpha */
        code = i << 4 * 3;
        StringInit ( & tok -> prefix, acc, 0, 0 );
        StringInit ( & tok -> alpha, acc, i, i );
    }

    /* remove alpha */
    acc += i;
    size -= i;

    /* scan digits */
    for ( i = 0; i < size; ++ i )
    {
        if ( ! isdigit ( acc [ i ] ) )
            break;
    }

    /* record digits */
    StringInit ( & tok -> digits, acc, i, i );
    StringInit ( & tok -> ext1, & acc [ i ], 0, 0 );
    tok -> ext2 = tok -> suffix = tok -> ext1;

    if ( i == 0 || i >= MAX_ACCESSION_LEN )
        return 0;

    code |= i << 4 * 2;

    /* normal return point for SRA */
    if ( i == size )
        return code;

    /* check for extension */
    if ( acc [ i ] != '.' )
        return 0;

    /* remove digit */
    acc += ++ i;
    size -= i;

    /* scan numeric extension */
    for ( i = 0; i < size; ++ i )
    {
        if ( ! isdigit ( acc [ i ] ) )
            break;
    }

    if ( i == 0 || i >= MAX_ACCESSION_LEN )
        return 0;

    /* record the actual extension */
    StringInit ( & tok -> ext1, acc, i, i );
    /* codify the extension simply as present, not by its length */
    code |= 1 << 4 * 1;

    if ( i == size )
        return code;

    /* scan for suffix */
    if ( acc [ i ] == '_' )
    {
        acc += ++ i;
        size -= i;
        for ( i = 0; i < size; ++ i )
        {
            if ( ! isalpha ( acc [ i ] ) )
                break;
        }

        /* this has to end the whole thing */
        if ( i == 0 || i != size )
            return 0;

        StringInit ( & tok -> suffix, acc, i, i );
        /* NB - not incorporating suffix into code right now */
        return code;
    }

    if ( acc [ i ] != '.' )
        return 0;


    /* remove first extension */
    acc += ++ i;
    size -= i;

    /* scan 2nd numeric extension */
    for ( i = 0; i < size; ++ i )
    {
        if ( ! isdigit ( acc [ i ] ) )
            break;
    }

    if ( i == 0 || i >= MAX_ACCESSION_LEN )
        return 0;

    StringInit ( & tok -> ext2, acc, i, i );
    code |= 1 << 4 * 0;

    if ( i == size )
        return code;

    /* no other patterns are accepted */
    return 0;
#endif
}


/* get_accession_app
 */
static
VResolverAppID get_accession_app ( const String * accession, bool refseq_ctx,
    VResolverAccToken *tok, bool *legacy_wgs_refseq )
{
    VResolverAppID app;
    uint32_t code = get_accession_code ( accession, tok );

    if (accession != NULL &&
        accession->addr != NULL && isdigit(accession->addr[0]))
    {
        /* TODO: KART */
        return appAny;
    }

    /* disregard extensions at this point */
    switch ( code >> 4 * 2 )
    {
    case 0x015: /* e.g. "J01415" or "J01415.2"     */
    case 0x026: /* e.g. "CM000071" or "CM000039.2" */
    case 0x126: /* e.g. "NZ_DS995509.1"            */
        app = appREFSEQ;
        break;

    case 0x036: /* e.g. "SRR012345"    */
    case 0x037: /* e.g. "SRR0123456"   */
    case 0x038: /* e.g. "SRR01234567"  */
    case 0x039: /* e.g. "SRR012345678" */

        /* detect accession with extension */
        if ( ( code & 0xFF ) != 0 )
            app = appAny;
        else
            app = appSRA;
        break;

    case 0x106: /* e.g. "NC_000012.10"                      */
    case 0x109: /* e.g. "NW_003315935.1", "GPC_000000393.1" */
        if ( tok -> prefix . size == 3 &&
             tok -> prefix . addr [ 0 ] == 'G' &&
             tok -> prefix . addr [ 1 ] == 'C' &&
             ( tok -> prefix . addr [ 2 ] == 'A' || tok -> prefix . addr [ 2 ] == 'F' ) )
        {
            /* e.g. "GCA_000392025.1_L" */
            app = appNAKMER;
            break;
        }

        app = appREFSEQ;
        break;

    case 0x042: /* e.g. "AAAB01" is WGS package name */
    case 0x048: /* e.g. "AAAA01000001"               */
    case 0x049: /* contig can be 6 or 7 digits       */
    case 0x142: /* e.g. "NZ_AAEW01"                  */
    case 0x148: /* e.g. "NZ_AAEW01000001"            */
    case 0x149:
        app = appWGS;
        break;

    case 0x029: /* e.g. NA000008777.1 */
        if ( code == 0x02910 )
        {
            if ( tok -> alpha . addr [ 0 ] == 'N' && tok -> alpha . addr [ 1 ] == 'A' )
            {
                app = appNANNOT;
                break;
            }
        }

        /* no break */

    default:
        /* TBD - people appear to be able to throw anything into refseq,
           so anything unrecognized we may as well test out there...
           but this should not stay the case */
        app = appREFSEQ;
    }

    if ( app == appWGS )
    {
        /* if we know this is for refseq, clobber it here */
        if ( refseq_ctx )
        {
            app = appREFSEQ;
            * legacy_wgs_refseq = true;
        }
    }

    return app;
}


/* LocalResolve
 *  resolve an accession into a VPath or not found
 *
 *  1. determine the type of accession we have, i.e. its "app"
 *  2. search all local algorithms of app type for accession
 *  3. return not found or new VPath
 */
static
rc_t VResolverLocalResolve ( const VResolver *self,
    const String * accession, const VPath ** path, bool refseq_ctx )
{
    rc_t rc;
    uint32_t i, count;

    VResolverAccToken tok;
    bool legacy_wgs_refseq = false;
    VResolverAppID app = get_accession_app ( accession, refseq_ctx, & tok, & legacy_wgs_refseq );

    /* search all local volumes by app and accession algorithm expansion */
    count = VectorLength ( & self -> local );
    for ( i = 0; i < count; ++ i )
    {
        const VResolverAlg *alg = VectorGet ( & self -> local, i );
        if ( alg -> app_id == app )
        {
            const bool for_cache = false;
            rc = VResolverAlgLocalResolve ( alg, self -> wd, & tok, path, legacy_wgs_refseq, for_cache );
            if ( rc == 0 )
                return 0;
        }
    }

    return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}

static
bool VPathHasRefseqContext ( const VPath * accession )
{
    size_t num_read;
    char option [ 64 ];
    rc_t rc = VPathOption ( accession, vpopt_vdb_ctx, option, sizeof option, & num_read );
    if ( rc != 0 )
        return false;
    return ( num_read == sizeof "refseq" - 1 &&
             strcase_cmp ( "refseq", sizeof "refseq" - 1,
                           option, num_read, num_read ) == 0 );
}


/* Local - DEPRECATED
 *  Find an existing local file/directory that is named by the accession.
 *  rcState of rcNotFound means it does not exist.
 *
 *  other rc code for failure are possible.
 *
 *  Accession must be an ncbi-acc scheme or a simple name with no 
 *  directory paths.
 */
LIB_EXPORT
rc_t CC VResolverLocal ( const VResolver * self,
    const VPath * accession, const VPath ** path )
{
    rc_t rc =  VResolverQuery ( self, eProtocolHttp, accession, path, NULL, NULL );
    if ( rc == 0 )
    {
        switch ( accession -> path_type )
        {
        case vpOID:
        case vpAccession:
        case vpNameOrOID:
        case vpNameOrAccession:
            if ( * path != accession )
                break;
        default:
            VPathRelease ( * path );
            * path = NULL;
            rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );
        }
    }
    return rc;
}


/* LocalEnable
 *  modify settings for using local repositories,
 *  meaning site, user-public and user-protected.
 *
 *  "enable" [ IN ] - enable or disable local access,
 *  or follow settings in KConfig
 *
 *  returns the previous state of "local-enabled" property
 *
 * NB - in VDB-2, the state is associated with library code
 *  shared libraries in separate closures will have separate
 *  state. this can only occur if dynamic ( manual ) loading of
 *  shared libraries is used, and will not occur with normal
 *  usage. in VDB-3 the state will be stored with the process,
 *  not the library.
 */
LIB_EXPORT
VResolverEnableState CC VResolverLocalEnable ( const VResolver * self, VResolverEnableState enable )
{
    int32_t val, cur, prior;

    if ( self == NULL )
        return vrUseConfig;

    /* convert "VResolverEnableState" to 32-bit signed integer for atomic operation */
    val = ( int32_t ) enable;

    /* before performing atomic swap, get the current setting,
       and return right away if it is already set correctly */
    prior = atomic32_read ( & enable_local );
    if ( prior != val ) do
    {
        cur = prior;
        prior = atomic32_test_and_set ( & enable_local, val, prior );
    }
    while ( prior != cur );

    return prior;
}


/* RemoteEnable
 *  apply or remove a process-wide enabling of remote access
 *  regardless of configuration settings
 *
 *  "enable" [ IN ] - if "true", enable all remote access
 *  if false, use settings from configuration.
 *
 *  returns the previous state of "remote-enabled" property
 *
 * NB - in VDB-2, the state is associated with library code
 *  shared libraries in separate closures will have separate
 *  state. this can only occur if dynamic ( manual ) loading of
 *  shared libraries is used, and will not occur with normal
 *  usage. in VDB-3 the state will be stored with the process,
 *  not the library.
 */
LIB_EXPORT
VResolverEnableState CC VResolverRemoteEnable ( const VResolver * self, VResolverEnableState enable )
{
    int32_t val, cur, prior;

    if ( self == NULL )
        return vrUseConfig;

    /* convert "VResolverEnableState" to 32-bit signed integer for atomic operation */
    val = ( int32_t ) enable;

    /* before performing atomic swap, get the current setting,
       and return right away if it is already set correctly */
    prior = atomic32_read ( & enable_remote );
    if ( prior != val ) do
    {
        cur = prior;
        prior = atomic32_test_and_set ( & enable_remote, val, prior );
    }
    while ( prior != cur );

    return prior;
}


/* CacheEnable
 *  modify settings for caching files in user repositories
 *
 *  "enable" [ IN ] - enable or disable user repository cache,
 *  or follow settings in KConfig
 *
 *  returns the previous state of "cache-enabled" property
 *
 * NB - in VDB-2, the state is associated with library code
 *  shared libraries in separate closures will have separate
 *  state. this can only occur if dynamic ( manual ) loading of
 *  shared libraries is used, and will not occur with normal
 *  usage. in VDB-3 the state will be stored with the process,
 *  not the library.
 */
LIB_EXPORT
VResolverEnableState CC VResolverCacheEnable ( const VResolver * self, VResolverEnableState enable )
{
    int32_t val, cur, prior;

    if ( self == NULL )
        return vrUseConfig;

    /* convert "VResolverEnableState" to 32-bit signed integer for atomic operation */
    val = ( int32_t ) enable;

    /* before performing atomic swap, get the current setting,
       and return right away if it is already set correctly */
    prior = atomic32_read ( & enable_cache );
    if ( prior != val ) do
    {
        cur = prior;
        prior = atomic32_test_and_set ( & enable_cache, val, prior );
    }
    while ( prior != cur );

    return prior;
}


/* RemoteResolve
 *  resolve an accession into a remote VPath or not found
 *  may optionally open a KFile to the object in the process
 *  of finding it
 *
 *  2. determine the type of accession we have, i.e. its "app"
 *  3. search all local algorithms of app type for accession
 *  4. return not found or new VPath
 */
static
rc_t VResolverRemoteResolve ( const VResolver *self,
    VRemoteProtocols protocols, const String * accession,
    const VPath ** path, const VPath **mapping,
    const KFile ** opt_file_rtn, bool refseq_ctx )
{
    rc_t rc, try_rc;
    uint32_t i, count;

    VResolverAccToken tok;
    VResolverAppID app, wildCard;
    bool legacy_wgs_refseq = false;

    VResolverEnableState remote_state = atomic32_read ( & enable_remote );

    /* subject the accession to pattern recognition */
    app = get_accession_app ( accession, refseq_ctx, & tok, & legacy_wgs_refseq );

    /* search all remote volumes by app and accession algorithm expansion */
    count = VectorLength ( & self -> remote );

    /* allow matching wildcard app */
    wildCard = appAny;
#if NO_REFSEQ_CGI
    if ( app == appREFSEQ )
        wildCard = -1;
#endif

    /* no error recorded yet */
    rc = 0;

    /* TBD - determine whether these settings interfere with
       case of resolving oid to cache location */

    /* test for forced enable, which applies only to main guys
       TBD - limit to main sub-category */
    if ( remote_state == vrAlwaysEnable )
    {
        for ( i = 0; i < count; ++ i )
        {
            const VResolverAlg *alg = VectorGet ( & self -> remote, i );
            if ( alg -> app_id == app || alg -> app_id == wildCard )
            {
                try_rc = VResolverAlgRemoteResolve ( alg, self -> kns, protocols, & tok, path, mapping, opt_file_rtn, legacy_wgs_refseq );
                if ( try_rc == 0 )
                    return 0;
                if ( rc == 0 )
                    rc = try_rc;
            }
        }
    }
    else
    {
        for ( i = 0; i < count; ++ i )
        {
            const VResolverAlg *alg = VectorGet ( & self -> remote, i );
            if ( ( alg -> app_id == app || alg -> app_id == wildCard ) && ! alg -> disabled )
            {
                try_rc = VResolverAlgRemoteResolve ( alg, self -> kns, protocols, & tok, path, mapping, opt_file_rtn, legacy_wgs_refseq );
                if ( try_rc == 0 )
                    return 0;
                if ( rc == 0 )
                    rc = try_rc;
            }
        }
    }

    if ( rc != 0 )
        return rc;

    return RC ( rcVFS, rcResolver, rcResolving, rcName, rcNotFound );
}


/* Remote
 *  Find an existing remote file that is named by the accession.
 *
 *  rcState of rcNotFound means it did not exist and can not be 
 *  downloaded. Probably a bad accession name.
 *
 *  Need a specific rc for no network configured.
 *  Need a specific rc for network access permitted.
 *
 *  Other rc code for failure are possible.
 *
 *  Accession must be an ncbi-acc scheme or a simple name with no 
 *  directory paths.
 *
 *  "opt_file_rtn" [ OUT, NULL OKAY ] - optional return parameter for
 *   any KFile that may be opened as a result of resolution. This can
 *   happen if resolving an accession involves opening a URL to a
 *   remote server, for example, in which case the KFile can be returned.
 */
LIB_EXPORT
rc_t CC VResolverRemote ( const VResolver * self,
    VRemoteProtocols protocols, const VPath * accession,
    const VPath ** path, const KFile ** opt_file_rtn )
{
    rc_t rc;

    if ( opt_file_rtn != NULL )
        * opt_file_rtn = NULL;

    rc = VResolverQuery ( self, protocols, accession, NULL, path, NULL );
    if ( rc == 0 && opt_file_rtn != NULL &&
        VPathGetUri_t ( * path ) != vpuri_fasp )
    {
#if USE_CURL
        char path_str [ 8192 ];
        rc = VPathReadUri ( * path, path_str, sizeof path_str, NULL );
        if ( rc == 0 )
            rc = KCurlFileMake ( opt_file_rtn, path_str, false);
        if ( rc != 0 )
        {
            VPathRelease ( * path );
            * path = NULL;
        }
#endif
    }

    return rc;
}


/* ExtractAccessionApp
 *  examine a path for accession portion,
 *  and try to recognize what app it belongs to
 */
static
VResolverAppID VResolverExtractAccessionApp ( const VResolver *self,
    const VPath * query, bool has_fragment,
    String * accession, VResolverAccToken * tok,
    bool *legacy_wgs_refseq )
{
    bool refseq_ctx = has_fragment;

    * accession = query -> path;

    if ( query -> fragment . size > 1 )
        refseq_ctx = true;

    /* should have something looking like an accession.
       determine its app to see if we were successful */
    return get_accession_app ( accession, refseq_ctx, tok, legacy_wgs_refseq );
}

static
bool VPathHasDownloadTicket ( const VPath * url )
{
    size_t num_read;
    char option [ 64 ];
    rc_t rc = VPathOption ( url, vpopt_gap_ticket, option, sizeof option, & num_read );
    return rc == 0;
}

static
rc_t VPathExtractAcc ( const VPath * url, VPath ** acc )
{
    rc_t rc;
    String accession;

    /* locate last path or accession guy */
    const char * start = string_rchr ( url -> path . addr, url -> path . size, '/' );
    const char * sep, * end = url -> path . addr + url -> path . size;
    if ( start ++ == NULL )
        start = url -> path . addr;

    /* strip off known extensions */
    sep = string_rchr ( start, end - start, '.' );
    while ( sep != NULL )
    {
        switch ( end - sep )
        {
        case 4:
            if ( strcase_cmp ( ".sra", 4, sep, 4, 4 ) == 0 )
                end = sep;
            else if ( strcase_cmp ( ".wgs", 4, sep, 4, 4 ) == 0 )
                end = sep;
            break;
        case 9:
            if ( strcase_cmp ( ".ncbi_enc", 9, sep, 9, 9 ) == 0 )
            {
                end = sep;
                sep = string_rchr ( start, end - start, '.' );
                continue;
            }
            break;
        }
        break;
    }

    /* this is the string */
    StringInit ( & accession, start, end - start, string_len ( start, end - start ) );

    /* now extract the mapping */
    rc = LegacyVPathMakeFmt ( acc, "ncbi-acc:%S%S%S",
        & accession, & url -> query, & url -> fragment );
    if ( rc == 0 )
    {
        VPath * ap = * acc;

        /* fix up case where we said accession but it was really a name */
        if ( ap -> acc_code == 0 || ap -> path_type != vpAccession )
            CONST_STRING ( & ap -> scheme, "ncbi-file" );
    }

    return rc;
}

static
rc_t VResolverCacheResolve ( const VResolver *self,
    const VPath * query, bool has_fragment,
    const VPath ** cache, bool refseq_ctx )
{
    rc_t rc = 0;

    String accession;
    VResolverAccToken tok;
    bool legacy_wgs_refseq = false;
    VResolverAppID app = VResolverExtractAccessionApp ( self,
        query, has_fragment, & accession, & tok, & legacy_wgs_refseq );

    /* going to walk the local volumes, and remember
       which one was best. actually, we have no algorithm
       for determining it, so it's just the comment for TBD */
    const VResolverAlg *alg, *best = NULL;

    /* search the volumes for a cache-enabled place */
    uint32_t i, count = VectorLength ( & self -> local );

    /* check for protected status by presence of a download ticket */
    bool protected = VPathHasDownloadTicket ( query );

    VResolverEnableState cache_state = atomic32_read ( & enable_cache );

    /* check for cache-enable override */
    if ( cache_state == vrAlwaysEnable )
    {
        for ( i = 0; i < count; ++ i )
        {
            alg = VectorGet ( & self -> local, i );
            if ( alg -> cache_capable && alg -> protected == protected &&
                 ( alg -> app_id == app || alg -> app_id == appAny ) )
            {
                /* try to find an existing cache file
                   NB - race condition exists unless
                   we do something with lock files */
                rc = VResolverAlgCacheResolve ( alg, self -> wd, & tok, cache, legacy_wgs_refseq );
                if ( rc == 0 )
                    return 0;

                /* just remember the first as best for now */
                if ( best == NULL )
                    best = alg;
            }
        }
    }
    else
    {
        for ( i = 0; i < count; ++ i )
        {
            alg = VectorGet ( & self -> local, i );
            if ( alg -> cache_enabled && alg -> protected == protected &&
                 ( alg -> app_id == app || alg -> app_id == appAny ) )
            {
                /* try to find an existing cache file
                   NB - race condition exists unless
                   we do something with lock files */
                rc = VResolverAlgCacheResolve ( alg, self -> wd, & tok, cache, legacy_wgs_refseq );
                if ( rc == 0 )
                    return 0;

                /* just remember the first as best for now */
                if ( best == NULL )
                    best = alg;
            }
        }
    }
    
    /* no existing cache file was found,
       so create a new one using the best
       TBD - this should remember a volume path */
    if ( best == NULL )
        rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );
    else
        rc = VResolverAlgMakeCachePath ( best, & tok, cache, legacy_wgs_refseq );

    return rc;
}


/* Cache
 *  Find a cache directory that might or might not contain a partially
 *  downloaded file.
 *
 *  Accession must be an ncbi-acc scheme, an http url or a simple name with no 
 *  directory paths. All three should return the same directory URL as a VPath. (?)
 *  Or should it be a directory or a file url depending upon finding a partial
 *  download? This would require co-ordination with all download mechanisms that
 *  we permit.
 *
 *  With refseq holding wgs objects we have a case were the downloaded file is not
 *  named the same as the original accession as the file archive you want is a
 *  container for other files.
 *
 *  Find local will give a path that has a special scheme in these cases. 
 *  Find remote will give the url for the container that contains the accession
 *  so using the returned VPath from resolve remote is better than the original
 *  accession in this one case.  I think...
 */
LIB_EXPORT
rc_t CC VResolverCache ( const VResolver * self,
    const VPath * url, const VPath ** path, uint64_t file_size )
{
    return VResolverQuery ( self, eProtocolHttp, url, NULL, NULL, path );
}

/* QueryOID
 */

static
rc_t get_query_accession ( const VPath * query, String * accession, char * oid_str, size_t bsize )
{
    rc_t rc;

    /* going to treat oid as accession */
    * accession = query -> path;

    /* if the VPath already gives us a numeral, great */
    if ( query -> path . size != 0 && query -> path . addr [ 0 ] != '0' )
        return 0;

    /* otherwise, generate one on stack */
    rc = string_printf ( oid_str, bsize, & accession -> size, "%u", query -> obj_id );
    if ( rc == 0 )
    {
        accession -> addr = oid_str;
        accession -> len = ( uint32_t ) accession -> size;
    }

    return rc;
}

static
rc_t VResolverQueryOID ( const VResolver * self, VRemoteProtocols protocols,
    const VPath * query, const VPath ** local, const VPath ** remote, const VPath ** cache )
{
    rc_t rc;

    /* require non-zero oid */
    if ( query -> obj_id == 0 )
        rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcCorrupt );
    else
    {
        /* temporary - no access to vfs
           NB - this manager will either use a singleton
           or create a new one with its existing config */
        VFSManager * vfs;
        rc = VFSManagerMake ( & vfs );
        if ( rc == 0 )
        {
            char oid_str [ 32 ];
            String accession;
            VPath * mapped_query = NULL;

            /* not expected to ever be true */
            bool refseq_ctx = VPathHasRefseqContext ( query );

            /* PREFACE - having an oid, we will need to map it to either
               an accession or simple filename before resolving to a
               local or cache path. there are two ways of getting this
               mapping: either through the VFS manager, or by asking the
               remote resolver CGI.

               ASSUMPTION - if the file exists locally or is cached,
               there should be a mapping available to VFS manager. this
               assumption can fail if the mapping database has been lost
               or damaged.
            */

            /* MAP OID TO ACCESSION */
            if ( local != NULL || cache != NULL )
            {
                /* we want a mapping. ask VFS manager for one */
                rc = VFSManagerGetObject ( vfs, query -> obj_id, & mapped_query );
                if ( GetRCState ( rc ) == rcNotFound )
                {
                    /* no mapping could be found. another possibility is to resolve remotely */
                    if ( remote != NULL || atomic32_read ( & enable_remote ) != vrAlwaysDisable )
                    {
                        rc = get_query_accession ( query, & accession, oid_str, sizeof oid_str );
                        if ( rc == 0 )
                        {
                            const VPath * remote2, * remote_mapping = NULL;
                            rc = VResolverRemoteResolve ( self, protocols, & accession,
                                & remote2, & remote_mapping, NULL, refseq_ctx );
                            if ( rc == 0 )
                            {
                                /* got it. now enter into VFS manager's table */
                                rc = VFSManagerRegisterObject ( vfs, query -> obj_id, remote_mapping );
                                if ( rc == 0 )
                                {
                                    mapped_query = ( VPath* ) remote_mapping;
                                    remote_mapping = NULL;
                                    if ( remote != NULL )
                                    {
                                        * remote = remote2;
                                        remote2 = NULL;
                                    }
                                }

                                VPathRelease ( remote2 );
                                VPathRelease ( remote_mapping );
                            }
                        }
                    }
                }

                if ( rc == 0 )
                {
                    assert ( mapped_query != NULL );

                    /* the returned VPath should be of a usable type */
                    assert ( mapped_query -> path_type == vpAccession       ||
                             mapped_query -> path_type == vpNameOrAccession ||
                             mapped_query -> path_type == vpName );
                    assert ( mapped_query -> path . size != 0 );
                }
            }

            /* RESOLVE FOR LOCAL PATH */
            if ( local != NULL && mapped_query != NULL )
            {
                /* grab the path as accession */
                accession = mapped_query -> path;

                /* resolve from accession to local path
                   will NOT find partial cache files */
                rc = VResolverLocalResolve ( self, & accession, local, refseq_ctx );
                if ( rc == 0 && remote != NULL && * remote != NULL )
                {
                    /* dump remote path used to map oid */
                    VPathRelease ( * remote );
                    * remote = NULL;
                }
            }

            if ( local == NULL || * local == NULL )
            {
                bool has_fragment = false;

                /* RESOLVE FOR REMOTE */
                if ( remote != NULL && * remote == NULL )
                {
                    rc = get_query_accession ( query, & accession, oid_str, sizeof oid_str );
                    if ( rc == 0 )
                    {
                        const VPath * remote_mapping = NULL;
                        rc = VResolverRemoteResolve ( self, protocols, & accession, remote,
                            ( mapped_query == NULL && cache != NULL ) ? & remote_mapping : NULL,
                            NULL, refseq_ctx );

                        if ( rc == 0 && mapped_query == NULL && cache != NULL && remote_mapping == NULL )
                        {
                            /* THIS IS LIKELY AN INTERNAL ERROR
                               EITHER THE CGI DID NOT RETURN A MAPPING
                               OR WE DID NOT PROPERLY PARSE IT */
                            VPathRelease ( * remote );
                            rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNull );
                        }

                        /* register new mapping */
                        if ( rc == 0 )
                        {
                            assert ( * remote != NULL );
                            if ( ( * remote ) -> fragment . size != 0 )
                                has_fragment = true;
                            if ( remote_mapping != NULL )
                            {
                                rc = VFSManagerRegisterObject ( vfs, query -> obj_id, remote_mapping );
                                if ( rc == 0 )
                                {
                                    mapped_query = ( VPath* ) remote_mapping;
                                    remote_mapping = NULL;
                                }
                                VPathRelease ( remote_mapping );
                            }
                        }
                    }
                }

                /* RESOLVE FOR CACHE */
                if ( ( remote == NULL || * remote != NULL ) && cache != NULL && mapped_query != NULL )
                {
                    /* resolve from accession to cache path */
                    rc = VResolverCacheResolve ( self, mapped_query, has_fragment, cache, refseq_ctx );
                    if ( rc != 0 && remote != NULL )
                    {
                        assert ( * cache == NULL );
                        VPathRelease ( * remote );
                        * remote = NULL;
                    }
                }
            }

            VPathRelease ( mapped_query );

            VFSManagerRelease ( vfs );
        }
    }

    return rc;
}

/* QueryAcc
 */
static
rc_t VResolverQueryAcc ( const VResolver * self, VRemoteProtocols protocols,
    const VPath * query, const VPath ** local, const VPath ** remote, const VPath ** cache )
{
    rc_t rc = 0;

    /* the accession should be directly usable */
    const String * accession = & query -> path;

    /* check if it is intended to locate a legacy refseq object */
    bool refseq_ctx = VPathHasRefseqContext ( query );

    /* will be needed to consult CGI */
    const VPath * remote2 = NULL, * mapped_query = NULL;

    /* LOCAL RESOLUTION */
    if ( local != NULL )
        rc = VResolverLocalResolve ( self, accession, local, refseq_ctx );

    if ( local == NULL || * local == NULL )
    {
        bool has_fragment = false;

        /* REMOTE RESOLUTION */
        if ( remote != NULL || ( self -> ticket != NULL && cache != NULL ) )
        {
            /* will need to map if protected */
            const VPath ** mapped_ptr = ( self -> ticket != NULL && cache != NULL ) ?
                & mapped_query : NULL;

            /* request remote resolution
               this does not need to map the query to an accession */
            rc = VResolverRemoteResolve ( self, protocols, accession,
                & remote2, mapped_ptr, NULL, refseq_ctx );

            if ( rc == 0 )
            {
                if ( remote2 -> fragment . size != 0 )
                    has_fragment = true;

                if ( remote != NULL )
                    * remote = remote2;
                else
                    VPathRelease ( remote2 );

                remote2 = NULL;
            }
        }

        if ( ( remote == NULL || * remote != NULL ) && cache != NULL )
        {
            if ( mapped_query != NULL )
                rc = VResolverCacheResolve ( self, mapped_query, has_fragment, cache, refseq_ctx );
#if 0
            /* the bad assumption that every remotely retrieved accession MUST be mapped */
            else if ( self -> ticket != NULL )
                rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );
#endif
            else
                rc = VResolverCacheResolve ( self, query, has_fragment, cache, refseq_ctx );

            if ( rc != 0 && remote != NULL )
            {
                assert ( * cache == NULL );
                if ( GetRCState ( rc ) == rcNotFound )
                    rc = 0;
                else
                {
                    VPathRelease ( * remote );
                    * remote = NULL;
                }
            }
        }

        if ( mapped_query != NULL )
            VPathRelease ( mapped_query );
    }

    return rc;
}

/* QueryPath
 *  this behavior may not be correct
 *  perhaps we should reject paths upon input,
 *  and only resolve things that need resolving
 *  but there is a thought that we can also transform paths
 */
static
rc_t VResolverQueryPath ( const VResolver * self, const VPath * query, const VPath ** local )
{
    rc_t rc;

    if ( local == NULL )
        return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );

    switch ( KDirectoryPathType ( self -> wd, "%.*s", ( int ) query -> path . size, query -> path . addr ) )
    {
    case kptFile:
    case kptDir:
    case kptCharDev:
    case kptBlockDev:
    case kptFIFO:
    case kptFile | kptAlias:
    case kptDir | kptAlias:
    case kptCharDev | kptAlias:
    case kptBlockDev | kptAlias:
    case kptFIFO | kptAlias:
        break;
    default:
        return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );
    }

    rc = VPathAddRef ( query );
    if ( rc == 0 )
        * local = query;

    return rc;
}


/* QueryName
 *  may eventually look for the name in local cache,
 *  but for now just return it as a path
 */
static
rc_t VResolverQueryName ( const VResolver * self, VRemoteProtocols protocols,
    const VPath * query, const VPath ** local, const VPath ** remote, const VPath ** cache )
{
    return VResolverQueryPath ( self, query, local );
}


/* QueryURL
 *  URL resolves to itself for remote and potentially to a path for cache
 */
static
rc_t VResolverQueryURL ( const VResolver * self, VRemoteProtocols protocols,
    const VPath * query, const VPath ** remote, const VPath ** cache )
{
    rc_t rc = 0;

    /* if neither remote nor cache, then must have requested local,
       and a URL cannot be resolved to local in our world... */
    if ( ( ( size_t ) remote | ( size_t ) cache ) == 0 )
        return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcIncorrect );

    /* the URL always resolves to itself for remote */
    if ( remote != NULL )
    {
        rc = VPathAddRef ( query );
        if ( rc != 0 )
            return rc;
        * remote = query;
    }

    /* if we want a cache location, then try to resolve it */
    if ( cache != NULL )
    {
        VPath *mapping;

        /* check for refseq context */
        bool refseq_ctx = VPathHasRefseqContext ( query );

        /* first, extract accession or name from URL */
        rc = VPathExtractAcc ( query, & mapping );
        if ( rc == 0 )
        {
            /* now map to cache location */
            rc = VResolverCacheResolve ( self, mapping, false, cache, refseq_ctx );
            VPathRelease ( mapping );
            if ( GetRCState ( rc ) == rcNotFound && remote != NULL )
                rc = 0;
        }

        /* any error must invalidate "remote" */
        if ( rc != 0 && remote != NULL )
        {
            VPathRelease ( * remote );
            * remote = NULL;
        }
    }

    return rc;
}


/* Query
 *  resolve object location to either an existing local path,
 *  or a pair of remote URL + local cache location.
 *
 *  "protocols" [ IN ] - the desired protocols for remote resolution
 *
 *  "query" [ IN ] - a path that can represent:
 *     accession : a recognizable accession from NCBI or known organization
 *     obj-id    : a dbGaP object id
 *     path      : a filesystem path
 *     url       : a remote location
 *
 *  "local" [ OUT, NULL OKAY ] - optional return parameter for local path:
 *     accession : resolve to local user or site path
 *     obj-id    : resolve to local user protected path
 *     path      : return duplicate of input
 *     url       : set to NULL
 *
 *  "remote" [ OUT, NULL OKAY ] - optional return parameter for remote path:
 *     accession : resolve to URL
 *     obj-id    : resolve to URL
 *     path      : set to NULL
 *     url       : set to duplicate
 *
 *  "cache" [ OUT, NULL OKAY ] - optional return parameter for cache path:
 *     accession : resolve to user cache path
 *     obj-id    : resolve to user cache path
 *     path      : set to NULL
 *     url       : resolve to user cache path
 *
 *  any of the output parameters may be NULL, but not all, i.e. there
 *  must be at least one non-NULL return parameter.
 *
 *  if you DON'T want local resolution, pass NULL for "local" and
 *  the query will be resolved remotely. if you don't want remote
 *  resolution, pass NULL for "remote".
 *
 *  a query that is resolved locally will always return NULL for
 *  "remote" and "cache", if the parameters are provided.
 */
LIB_EXPORT
rc_t CC VResolverQuery ( const VResolver * self, VRemoteProtocols protocols,
    const VPath * query, const VPath ** local, const VPath ** remote, const VPath ** cache )
{
    rc_t rc;

    if ( ( ( size_t ) local | ( size_t ) remote | ( size_t ) cache ) == 0 )
        rc = RC ( rcVFS, rcResolver, rcResolving, rcParam, rcNull );
    else
    {
        if ( local != NULL )
        {
            * local = NULL;
            if ( atomic32_read ( & enable_local ) == vrAlwaysDisable )
                local = NULL;
        }

        if ( remote != NULL )
        {
            * remote = NULL;
            if ( protocols >= eProtocolLastDefined )
                return RC ( rcVFS, rcResolver, rcResolving, rcParam, rcInvalid );
            if ( atomic32_read ( & enable_remote ) == vrAlwaysDisable )
                remote = NULL;
        }

        if ( cache != NULL )
        {
            * cache = NULL;
            if ( atomic32_read ( & enable_cache ) == vrAlwaysDisable )
                cache = NULL;
        }

        if ( self == NULL )
            rc = RC ( rcVFS, rcResolver, rcResolving, rcSelf, rcNull );
        else if ( query == NULL )
            rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNull );
        else if ( ( ( size_t ) local | ( size_t ) remote | ( size_t ) cache ) == 0 )
            rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcNotFound );
        else
        {
            switch ( query -> scheme_type )
            {
            case vpuri_none:
            case vpuri_ncbi_file:
            case vpuri_file:
            case vpuri_ncbi_acc:
            case vpuri_ncbi_obj:
                break;

            case vpuri_http:
                switch ( protocols )
                {
                case eProtocolHttp:
                case eProtocolFaspHttp:
                case eProtocolHttpFasp:
                    return VResolverQueryURL ( self, protocols, query, remote, cache );
                }
                return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcIncorrect );

            case vpuri_fasp:
                switch ( protocols )
                {
                case eProtocolFasp:
                case eProtocolFaspHttp:
                case eProtocolHttpFasp:
                    return VResolverQueryURL ( self, protocols, query, remote, cache );
                }
                return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcIncorrect );

            default:
                return RC ( rcVFS, rcResolver, rcResolving, rcPath, rcIncorrect );
            }

            switch ( query -> path_type )
            {
            case vpInvalid:
                rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcInvalid );
                break;

            case vpOID:
                rc = VResolverQueryOID ( self, protocols, query, local, remote, cache );
                break;

            case vpAccession:
                rc = VResolverQueryAcc ( self, protocols, query, local, remote, cache );
                break;

            case vpNameOrOID:
                rc = VResolverQueryOID ( self, protocols, query, local, remote, cache );
                if ( rc != 0 )
                    goto try_name;
                break;

            case vpNameOrAccession:
                rc = VResolverQueryAcc ( self, protocols, query, local, remote, cache );
                if ( rc != 0 )
                    goto try_name;
                break;

            case vpName:
            try_name:
                rc = VResolverQueryName ( self, protocols, query, local, remote, cache );
                break;

            case vpRelPath:
            case vpFullPath:
            case vpUNCPath:
                rc = VResolverQueryPath ( self, query, local );
                break;

            default:
                rc = RC ( rcVFS, rcResolver, rcResolving, rcPath, rcIncorrect );
            }
        }
    }

    return rc;
}



/* LoadVolume
 *  capture volume path and other information
 */
static
rc_t VResolverAlgLoadVolume ( VResolverAlg *self, uint32_t *num_vols, const char *start, size_t size )
{
    rc_t rc = 0;

#if 0
    /* trim volume whitespace */
    while ( size != 0 && isspace ( start [ 0 ] ) )
    {
        ++ start;
        -- size;
    }
    while ( size != 0 && isspace ( start [ size - 1 ] ) )
        -- size;
#endif

    /* trim trailing slashes */
    while ( size != 0 && start [ size - 1 ] == '/' )
        -- size;

    /* now see if the string survives */
    if ( size != 0 )
    {
        String loc_vol_str;
        const String *vol_str;
        StringInit ( & loc_vol_str, start, size, string_len ( start, size ) );
        rc = StringCopy ( & vol_str, & loc_vol_str );
        if ( rc == 0 )
        {
            rc = VectorAppend ( & self -> vols, NULL, vol_str );
            if ( rc == 0 )
            {
                * num_vols += 1;
                return 0;
            }

            StringWhack ( vol_str );
        }
    }

    return rc;
}

/* LoadVolumes
 *
 *    path-list
 *        = PATH
 *        | <path-list> ':' PATH ;
 */
static
rc_t VResolverAlgLoadVolumes ( VResolverAlg *self, uint32_t *num_vols, const String *vol_list )
{
    const char *start = vol_list -> addr;
    const char *end = & vol_list -> addr [ vol_list -> size ];
    const char *sep = string_chr ( start, end - start, ':' );
    while ( sep != NULL )
    {
        rc_t rc = VResolverAlgLoadVolume ( self, num_vols, start, sep - start );
        if ( rc != 0 )
            return rc;
        start = sep + 1;
        sep = string_chr ( start, end - start, ':' );
    }
    return VResolverAlgLoadVolume ( self, num_vols, start, end - start );
}

/* LoadAlgVolumes
 *
 *    volumes
 *        = <path-list> ;
 */
static
rc_t VResolverLoadAlgVolumes ( Vector *algs, const String *root, const String *ticket,
    VResolverCacheAllow allow_cache, VResolverAppID app_id, VResolverAlgID alg_id,
     uint32_t *num_vols, const String *vol_list, bool protected, bool disabled, bool caching )
{
    VResolverAlg *alg;
    rc_t rc = VResolverAlgMake ( & alg, root, app_id, alg_id, protected, disabled );
    if ( rc == 0 )
    {
        alg -> ticket = ticket;
        alg -> cache_capable = allow_cache == cacheAllow;
        alg -> cache_enabled = caching;

        if ( ticket != NULL )
            alg -> alg_id = algCGI;

        rc = VResolverAlgLoadVolumes ( alg, num_vols, vol_list );
        if ( rc == 0 && VectorLength ( & alg -> vols ) != 0 )
        {
            rc = VectorAppend ( algs, NULL, alg );
            if ( rc == 0 )
                return 0;
        }

        VResolverAlgWhack ( alg, NULL );
    }

    return rc;
}

/* LoadApp
 *
 *    alg-block
 *        = <alg-type> <volumes> ;
 *
 *    alg-type
 *        = "flat" | "sraFlat" | "sra1024" | "sra1000" | "fuse1000"
 *        | "refseq" | "wgs" | "wgsFlag" | "fuseWGS"
 *        | "ncbi" | "ddbj" | "ebi"
 *        | "nannot" | "nannotFlat" | "fuseNANNOT" ;
 */
static
rc_t VResolverLoadVolumes ( Vector *algs, const String *root, const String *ticket,
    VResolverCacheAllow allow_cache, VResolverAppID app_id, uint32_t *num_vols,
    const KConfigNode *vols, bool resolver_cgi, bool protected, bool disabled, bool caching )
{
    KNamelist *algnames;
    rc_t rc = KConfigNodeListChildren ( vols, & algnames );
    if ( rc == 0 )
    {
        uint32_t i, count;
        rc = KNamelistCount ( algnames, & count );
        for ( i = 0; i < count && rc == 0; ++ i )
        {
            const char *algname;
            rc = KNamelistGet ( algnames, i, & algname );
            if ( rc == 0 )
            {
                const KConfigNode *alg;
                rc = KConfigNodeOpenNodeRead ( vols, & alg, algname );
                if ( rc == 0 )
                {
                    VResolverAlgID alg_id = algUnknown;

                    /* if using CGI for resolution */
                    if ( resolver_cgi || strcmp ( algname, "cgi" ) == 0 )
                        alg_id = algCGI;
                    /* stored in a flat directory with ".sra" extension */
                    else if ( strcmp ( algname, "sraFlat" ) == 0 )
                        alg_id = algSRAFlat;
                    /* stored in a three-level directory with 1K banks and ".sra" extension */
                    else if ( strcmp ( algname, "sra1024" ) == 0 )
                        alg_id = algSRA1024;
                    /* stored in a three-level directory with 1000 banks and ".sra" extension */
                    else if ( strcmp ( algname, "sra1000" ) == 0 )
                        alg_id = algSRA1000;
                    /* stored in a four-level directory with 1000 banks and ".sra" extension */
                    else if ( strcmp ( algname, "fuse1000" ) == 0 )
                        alg_id = algFUSE1000;
                    /* stored in a flat directory with no extension */
                    else if ( strcmp ( algname, "refseq" ) == 0 )
                        alg_id = algREFSEQ;
                    /* stored in a flat directory with no extension */
                    else if ( strcmp ( algname, "wgsFlat" ) == 0 )
                        alg_id = algWGSFlat;
                    /* stored in a multi-level directory with no extension */
                    else if ( strcmp ( algname, "wgs" ) == 0 )
                        alg_id = algWGS;
                    else if ( strcmp ( algname, "fuseWGS" ) == 0 )
                        alg_id = algFuseWGS;
                    /* stored in a three-level directory with 1K banks and no extension */
                    else if ( strcmp ( algname, "ncbi" ) == 0 ||
                              strcmp ( algname, "ddbj" ) == 0 )
                        alg_id = algSRA_NCBI;
                    /* stored in a three-level directory with 1000 banks and no extension */
                    else if ( strcmp ( algname, "ebi" ) == 0 )
                        alg_id = algSRA_EBI;

                    /* new named annotation */
                    else if ( strcmp ( algname, "nannotFlat" ) == 0 )
                        alg_id = algNANNOTFlat;
                    else if ( strcmp ( algname, "nannot" ) == 0 )
                        alg_id = algNANNOT;
                    else if ( strcmp ( algname, "fuseNANNOT" ) == 0 )
                        alg_id = algFuseNANNOT;

                    /* new named annotation */
                    else if ( strcmp ( algname, "nakmerFlat" ) == 0 )
                        alg_id = algNAKMERFlat;
                    else if ( strcmp ( algname, "nakmer" ) == 0 )
                        alg_id = algNAKMER;
                    else if ( strcmp ( algname, "fuseNAKMER" ) == 0 )
                        alg_id = algFuseNAKMER;

                    if ( alg_id != algUnknown )
                    {
                        String *vol_list;
                        rc = KConfigNodeReadString ( alg, & vol_list );
                        if ( rc == 0 )
                        {
                            if ( StringLength ( vol_list ) != 0 )
                            {
                                rc = VResolverLoadAlgVolumes ( algs, root, ticket, allow_cache,
                                    app_id, alg_id, num_vols, vol_list, protected, disabled, caching );
                            }
                            StringWhack ( vol_list );
                        }
                    }

                    KConfigNodeRelease ( alg );
                }
            }
        }

        KNamelistRelease ( algnames );
    }
    return rc;
}

/* LoadApp
 *
 *    app
 *        = [ <disabled> ] [ <cache-enabled> ] <vol-group> ;
 *
 *    disabled
 *        = "disabled" '=' ( "true" | "false" ) ;
 *
 *    cache-enabled
 *        = "cache-enabled" '=' ( "true" | "false" ) ;
 *
 *    vol-group
 *        = "volumes" <alg-block>* ;
 */
static
rc_t VResolverLoadApp ( VResolver *self, Vector *algs, const String *root, const String *ticket,
    VResolverCacheAllow allow_cache, VResolverAppID app_id, uint32_t *num_vols,
    const KConfigNode *app, bool resolver_cgi, bool protected, bool disabled, bool caching )
{
    const KConfigNode *node;

    /* test for disabled app - it is entirely possible */
    rc_t rc = KConfigNodeOpenNodeRead ( app, & node, "disabled" );
    if ( rc == 0 )
    {
        bool app_disabled = false;
        rc = KConfigNodeReadBool ( node, & app_disabled );
        KConfigNodeRelease ( node );
        if ( rc == 0 && app_disabled && algs == & self -> local )
            return 0;
        disabled |= app_disabled;
    }

    /* test again for cache enabled */
    if ( allow_cache == cacheAllow )
    {
        rc = KConfigNodeOpenNodeRead ( app, & node, "cache-enabled" );
        if ( rc == 0 )
        {
            /* allow this node to override current value */
            bool cache;
            rc = KConfigNodeReadBool ( node, & cache );
            KConfigNodeRelease ( node );
            if ( rc == 0 )
                caching = cache;
        }
    }

    /* get volumes */
    rc = KConfigNodeOpenNodeRead ( app, & node, "volumes" );
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        rc = VResolverLoadVolumes ( algs, root, ticket, allow_cache,
            app_id, num_vols, node, resolver_cgi, protected, disabled, caching );
        KConfigNodeRelease ( node );
    }

    return rc;
}

/* LoadApps
 *
 *    app-block
 *        = <app-name> <app> ;
 *
 *    app-name
 *        = "refseq" | "sra" | "wgs" | "nannot" | "nakmer" ;
 */
static
rc_t VResolverLoadApps ( VResolver *self, Vector *algs, const String *root,
    const String *ticket, VResolverCacheAllow allow_cache, const KConfigNode *apps,
    bool resolver_cgi, bool protected, bool disabled, bool caching )
{
    KNamelist *appnames;
    rc_t rc = KConfigNodeListChildren ( apps, & appnames );
    if ( rc == 0 )
    {
        uint32_t i, count;
        rc = KNamelistCount ( appnames, & count );
        if ( resolver_cgi && rc == 0 && count == 0 )
        {
            VResolverAlg *cgi;
            rc = VResolverAlgMake ( & cgi, root, appAny, algCGI, protected, disabled );
            if ( rc == 0 )
            {
                rc = VectorAppend ( algs, NULL, cgi );
                if ( rc == 0 )
                {
                    ++ self -> num_app_vols [ appAny ];
                    return 0;
                }
            }
        }
        else for ( i = 0; i < count && rc == 0; ++ i )
        {
            const char *appname;
            rc = KNamelistGet ( appnames, i, & appname );
            if ( rc == 0 )
            {
                const KConfigNode *app;
                rc = KConfigNodeOpenNodeRead ( apps, & app, appname );
                if ( rc == 0 )
                {
                    VResolverAppID app_id = appUnknown;
                    if ( strcmp ( appname, "refseq" ) == 0 )
                        app_id = appREFSEQ;
                    else if ( strcmp ( appname, "sra" ) == 0 )
                        app_id = appSRA;
                    else if ( strcmp ( appname, "wgs" ) == 0 )
                        app_id = appWGS;
                    else if ( strcmp ( appname, "nannot" ) == 0 )
                        app_id = appNANNOT;
                    else if ( strcmp ( appname, "nakmer" ) == 0 )
                        app_id = appNAKMER;

                    rc = VResolverLoadApp ( self, algs, root, ticket, allow_cache, app_id,
                        & self -> num_app_vols [ app_id ], app, resolver_cgi, protected, disabled, caching );

                    KConfigNodeRelease ( app );
                }
            }
        }
        KNamelistRelease ( appnames );
    }
    return rc;
}

/* LoadRepo
 *
 *    repository
 *        = [ <disabled> ] [ <cache-enabled> ] <root> <app-group> ;
 *
 *    disabled
 *        = "disabled" '=' ( "true" | "false" ) ;
 *
 *    cache-enabled
 *        = "cache-enabled" '=' ( "true" | "false" ) ;
 *
 *    root
 *        = "root" '=' PATH ;
 *
 *    app-group
 *        = "apps" <app-block>* ;
 */
static
rc_t VResolverLoadRepo ( VResolver *self, Vector *algs, const KConfigNode *repo,
    const String *ticket, VResolverCacheAllow allow_cache, bool protected )
{
    const KConfigNode *node;
    bool caching, resolver_cgi;

    /* test for disabled repository */
    bool disabled = false;
    rc_t rc = KConfigNodeOpenNodeRead ( repo, & node, "disabled" );
    if ( rc == 0 )
    {
        rc = KConfigNodeReadBool ( node, & disabled );
        KConfigNodeRelease ( node );

        /* don't bother recording local, disabled repositories */
        if ( rc == 0 && disabled && algs == & self -> local )
            return 0;
    }

    /* check for caching */
    caching = allow_cache;
    if ( allow_cache )
    {
        rc = KConfigNodeOpenNodeRead ( repo, & node, "cache-enabled" );
        if ( rc == 0 )
        {
            rc = KConfigNodeReadBool ( node, & caching );
            KConfigNodeRelease ( node );
            if ( rc != 0 )
                caching = false;
        }
    }

    /* cache-capable repositories cannot be remote resolvers */
    resolver_cgi = false;
    if ( allow_cache )
        rc = KConfigNodeOpenNodeRead ( repo, & node, "root" );
    else
    {
        /* check for specific resolver CGI */
        rc = KConfigNodeOpenNodeRead ( repo, & node, "resolver-cgi" );
        if ( rc == 0 )
            resolver_cgi = true;
        /* or get the repository root */
        else if ( GetRCState ( rc ) == rcNotFound )
        {
            rc = KConfigNodeOpenNodeRead ( repo, & node, "root" );
        }
    }
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        /* read root as String */
        String *root;
        rc = KConfigNodeReadString ( node, & root );
        KConfigNodeRelease ( node );
        if ( GetRCState ( rc ) == rcNotFound )
            rc = 0;
        else if ( rc == 0 )
        {
            /* perform a bit of cleanup on root */
            while ( root -> size != 0 && root -> addr [ root -> size - 1 ] == '/' )
            {
                /* this is terribly nasty, but known to be safe */
                -- root -> len;
                -- root -> size;
            }

            /* store string on VResolver for management purposes,
               pass the loaned reference to sub-structures */
            rc = VectorAppend ( & self -> roots, NULL, root );
            if ( rc != 0 )
                StringWhack ( root );
            else
            {
                /* open the "apps" sub-node */
                rc = KConfigNodeOpenNodeRead ( repo, & node, "apps" );
                if ( rc == 0 )
                {
                    rc = VResolverLoadApps ( self, algs, root, ticket,
                        allow_cache, node, resolver_cgi, protected, disabled, caching );
                    KConfigNodeRelease ( node );
                }
                else if ( GetRCState ( rc ) == rcNotFound )
                {
                    rc = 0;
                    if ( resolver_cgi )
                    {
                        VResolverAlg *cgi;
                        rc = VResolverAlgMake ( & cgi, root, appAny, algCGI, protected, disabled );
                        if ( rc == 0 )
                        {
                            cgi -> ticket = ticket;

                            rc = VectorAppend ( algs, NULL, cgi );
                            if ( rc == 0 )
                            {
                                ++ self -> num_app_vols [ appAny ];
                                return 0;
                            }
                        }

                        VResolverAlgWhack ( cgi, NULL );
                    }
                }
            }
        }
    }

    return rc;
}


/* LoadNamedRepo
 *
 *    repository-block
 *        = ID <repository> ;
 */
static
rc_t VResolverLoadNamedRepo ( VResolver *self, Vector *algs, const KConfigNode *sub,
    const String *ticket, const char *name, VResolverCacheAllow allow_cache, bool protected )
{
    const KConfigNode *repo;
    rc_t rc = KConfigNodeOpenNodeRead ( sub, & repo, name );
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        rc = VResolverLoadRepo ( self, algs, repo, ticket, allow_cache, protected );
        KConfigNodeRelease ( repo );
    }
    return rc;
}

/* LoadSubCategory
 *
 *    sub-category-block
 *        = <sub-category> <repository-block>* ;
 *
 *    sub-category
 *        = "main" | "aux" | "protected"
 *
 *    repository-block
 *        = ID <repository> ;
 */
static
rc_t VResolverLoadSubCategory ( VResolver *self, Vector *algs, const KConfigNode *kfg,
    const String *ticket, const char *sub_path, VResolverCacheAllow allow_cache, bool protected )
{
    const KConfigNode *sub;
    rc_t rc = KConfigNodeOpenNodeRead ( kfg, & sub, sub_path );
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        KNamelist *children;
        rc = KConfigNodeListChildren ( sub, & children );
        if ( rc == 0 )
        {
            uint32_t i, count;
            rc = KNamelistCount ( children, & count );
            for ( i = 0; i < count && rc == 0; ++ i )
            {
                const char *name;
                rc = KNamelistGet ( children, i, & name );
                if ( rc == 0 )
                    rc = VResolverLoadNamedRepo ( self, algs, sub, ticket, name, allow_cache, protected );
            }

            KNamelistRelease ( children );
        }
        KConfigNodeRelease ( sub );
    }
    return rc;
}

/* LoadProtected
 *  special function to handle single, active protected workspace
 */
static
rc_t VResolverLoadProtected ( VResolver *self, const KConfigNode *kfg, const char *rep_name )
{
    const KConfigNode *repo;
    rc_t rc = KConfigNodeOpenNodeRead ( kfg, & repo, "user/protected/%s", rep_name );
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        rc = VResolverLoadRepo ( self, & self -> local, repo, NULL, cacheAllow, true );
        KConfigNodeRelease ( repo );
    }
    return rc;
}

/* LoadLegacyRefseq
 *  load refseq information from KConfig
 *
 *  there are two legacy versions being supported
 *
 *    legacy-refseq
 *        = "refseq" <legacy-vol-or-repo> ;
 *
 *    legacy-vol-or-repo
 *        = "volumes" '=' <path-list>
 *        | <legacy-refseq-repo> <legacy-refseq-vols>
 *        ;
 */
static
rc_t VResolverLoadLegacyRefseq ( VResolver *self, const KConfig *cfg )
{
    const KConfigNode *vols;
    rc_t rc = KConfigOpenNodeRead ( cfg, & vols, "/refseq/paths" );
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        String *vol_list;
        rc = KConfigNodeReadString ( vols, & vol_list );
        if ( rc == 0 )
        {
            const bool protected = false;
            const bool disabled = false;
            const bool caching = true;
            rc = VResolverLoadAlgVolumes ( & self -> local, NULL, NULL, cacheAllow,
                appREFSEQ, algREFSEQ,  & self -> num_app_vols [ appREFSEQ ],
                vol_list, protected, disabled, caching );
            StringWhack ( vol_list );
        }
        KConfigNodeRelease ( vols );
    }

    return rc;
}


/* ForceRemoteRefseq
 *  makes sure there is a remote source of refseq
 *  or else adds a hard-coded URL to NCBI
 */
static
rc_t VResolverForceRemoteRefseq ( VResolver *self )
{
    rc_t rc;
    bool found;
    String local_root;
    const String *root;

    uint32_t i, count = VectorLength ( & self -> remote );
    for ( found = false, i = 0; i < count; ++ i )
    {
        VResolverAlg *alg = ( VResolverAlg* ) VectorGet ( & self -> remote, i );
        if ( alg -> app_id == appREFSEQ )
        {
            found = true;
            if ( alg -> disabled )
                alg -> disabled = false;
        }
    }

    if ( found )
        return 0;

    if ( self -> num_app_vols [ appAny ] != 0 )
    {
        for ( i = 0; i < count; ++ i )
        {
            VResolverAlg *alg = ( VResolverAlg* ) VectorGet ( & self -> remote, i );
            if ( alg -> app_id == appAny )
            {
                found = true;
                if ( alg -> disabled )
                    alg -> disabled = false;
            }
        }
    }

    if ( found )
        return 0;

    /* create one from hard-coded constants */
    StringInitCString ( & local_root, "http://ftp-trace.ncbi.nlm.nih.gov/sra" );
    rc = StringCopy ( & root, & local_root );    
    if ( rc == 0 )
    {
        rc = VectorAppend ( & self -> roots, NULL, root );
        if ( rc != 0 )
            StringWhack ( root );
        else
        {
            String vol_list;
            const bool protected = false;
            const bool disabled = false;
            const bool caching = false;
            StringInitCString ( & vol_list, "refseq" );
            rc = VResolverLoadAlgVolumes ( & self -> remote, root, NULL, cacheDisallow,
                appREFSEQ, algREFSEQ, & self -> num_app_vols [ appREFSEQ ],
                & vol_list, protected, disabled, caching );
        }
    }

    return rc;
}


/* GetDownloadTicket
 *  if we are within a working environment that has a download ticket,
 *  capture it here and add that local repository into the mix
 */
static
const String *VResolverGetDownloadTicket ( const VResolver *self,
    const KRepository *protected, char *buffer, size_t bsize )
{
    const String *ticket = NULL;
    if ( protected != NULL )
    {
        rc_t rc = KRepositoryName ( protected, buffer, bsize, NULL );
        if ( rc == 0 )
        {
            size_t ticsz;
            char ticbuf [ 256 ];
            rc = KRepositoryDownloadTicket ( protected, ticbuf, sizeof ticbuf, & ticsz );
            if ( rc == 0 )
            {
                String tic;
                StringInit ( & tic, ticbuf, ticsz, ( uint32_t ) ticsz );
                rc = StringCopy ( & ticket, & tic );
            }
        }
    }
    return ticket;
}


/* ForceRemoteProtected
 *  makes sure there is a remote CGI
 */
static
rc_t VResolverForceRemoteProtected ( VResolver *self )
{
    rc_t rc;
    const String *root;

    /* create one from hard-coded constants */
    String cgi_root;
    StringInitCString ( & cgi_root, "http://www.ncbi.nlm.nih.gov/Traces/names/names.cgi" );
    rc = StringCopy ( & root, & cgi_root );    
    if ( rc == 0 )
    {
        rc = VectorAppend ( & self -> roots, NULL, root );
        if ( rc != 0 )
            StringWhack ( root );
        else
        {
            const bool protected = true;
            const bool disabled = false;

            VResolverAlg *cgi;
            rc = VResolverAlgMake ( & cgi, root, appAny, algCGI, protected, disabled );
            if ( rc == 0 )
            {
                cgi -> ticket = self -> ticket;

                rc = VectorAppend ( & self -> remote, NULL, cgi );
                if ( rc == 0 )
                {
                    ++ self -> num_app_vols [ appAny ];
                    return 0;
                }
            }

            VResolverAlgWhack ( cgi, NULL );
        }
    }

    return rc;
}


/* Load
 *  load the respository from ( current ) KConfig
 *
 *  using pseudo BNF, it looks like this:
 *
 *    repositories
 *        = "repository" <category-block>* ;
 *
 *    category-block
 *        = <category> <sub-category-block>* ;
 *
 *    category
 *        = "remote" | "site" | "user" ;
 *
 *    sub-category-block
 *        = <sub-category> <repository-block>* ;
 *
 *    sub-category
 *        = "main" | "aux" | "protected"
 */
static
rc_t VResolverDetectSRALeafPath ( VResolver *self )
{
    /* capture working directory as "root" path */
    const KDirectory *wd = self -> wd;
    char cwd [ 4096 ];
    rc_t rc = KDirectoryResolvePath ( wd, true, cwd, sizeof cwd, "." );
    if ( rc == 0 )
    {
        const String *root;

        /* convert C-string to real string */
        String cwd_str;
        StringInitCString ( & cwd_str, cwd );

        /* create a copy on heap */
        rc = StringCopy ( & root, & cwd_str );
        if ( rc == 0 )
        {
            /* insert into "roots" */
            rc = VectorAppend ( & self -> roots, NULL, root );
            if ( rc == 0 )
            {
                /* create an algorithm for any application where the
                   spec is to be treated as a leaf path */
                VResolverAlg *alg;
                rc = VResolverAlgMake ( & alg, root, appAny, algLeafPath, self -> ticket != NULL, false );
                if ( rc == 0 )
                {
                    const String *vol;

                    /* create a single volume - "." */
                    CONST_STRING ( & cwd_str, "." );
                    rc = StringCopy ( & vol, & cwd_str );
                    if ( rc == 0 )
                    {
                        rc = VectorAppend ( & alg -> vols, NULL, vol );
                        if ( rc != 0 )
                            free ( ( void* ) vol );
                        else
                        {
                            /* insert into local resolution path */
                            rc = VectorAppend ( & self -> local, NULL, alg );
                            if ( rc == 0 )
                                return 0;
                        }
                    }
                
                    VResolverAlgWhack ( alg, NULL );
                }
            }

            free ( ( void* ) root );
        }
    }
    return rc;
}

static
rc_t VResolverLoad ( VResolver *self, const KRepository *protected, const KConfig *cfg )
{
    bool have_remote_protected = false;

    const KConfigNode *kfg;
    rc_t rc = KConfigOpenNodeRead ( cfg, & kfg, "repository" );
    if ( GetRCState ( rc ) == rcNotFound )
        rc = 0;
    else if ( rc == 0 )
    {
        /* check to see what the current directory is */
        char buffer [ 256 ];
        self -> ticket = VResolverGetDownloadTicket ( self, protected, buffer, sizeof buffer );

        /* allow user to specify leaf paths in current directory */
        rc = VResolverDetectSRALeafPath ( self );

        /* if the user is inside of a protected workspace, load it now */
        if ( rc == 0 && self -> ticket != NULL )
            rc = VResolverLoadProtected ( self, kfg, buffer );

        /* now load user public repositories */
        if ( rc == 0 )
            rc = VResolverLoadSubCategory ( self, & self -> local, kfg, NULL, "user/main", cacheAllow, false );
        if ( rc == 0 )
            rc = VResolverLoadSubCategory ( self, & self -> local, kfg, NULL, "user/aux", cacheAllow, false );

        /* load any site repositories */
        if ( rc == 0 )
            rc = VResolverLoadSubCategory ( self, & self -> local, kfg, NULL, "site/main", cacheDisallow, false );
        if ( rc == 0 )
            rc = VResolverLoadSubCategory ( self, & self -> local, kfg, NULL, "site/aux", cacheDisallow, false );

        /* if within a protected workspace, load protected remote repositories */
        if ( rc == 0 && self -> ticket != NULL )
        {
            rc = KNSManagerMake ( ( KNSManager** ) & self -> kns );
            if ( rc == 0 )
            {
                uint32_t entry_vols = VectorLength ( & self -> remote );
                rc = VResolverLoadSubCategory ( self, & self -> remote, kfg,
                    self -> ticket, "remote/protected", cacheDisallow, true );
                have_remote_protected = VectorLength ( & self -> remote ) > entry_vols;
            }
        }

        /* load any remote repositories */
        if ( rc == 0 )
            rc = VResolverLoadSubCategory ( self, & self -> remote, kfg, NULL, "remote/main", cacheDisallow, false );
        if ( rc == 0 )
            rc = VResolverLoadSubCategory ( self, & self -> remote, kfg, NULL, "remote/aux", cacheDisallow, false );

        KConfigNodeRelease ( kfg );

        /* recover from public remote repositories using resolver CGI */
        if ( self -> kns == NULL
#if USE_CURL
             && self -> num_app_vols [ appAny ] != 0
#endif
            )
        {
            rc = KNSManagerMake ( ( KNSManager** ) & self -> kns );
        }
    }

    if ( rc == 0 && self -> num_app_vols [ appAny ] == 0 )
    {
        bool has_current_refseq = true;

        /* AT THIS POINT, a current configuration will have something.
           But, older out-of-date configurations may exist and need special handling. */
        if ( self -> num_app_vols [ appREFSEQ ] == 0 )
        {
            has_current_refseq = false;
            rc = VResolverLoadLegacyRefseq ( self, cfg );
        }

        /* now, one more special case - for external users
           who had legacy refseq configuration but nothing for SRA,
           force existence of a remote refseq access */
        if ( rc == 0
             && ! has_current_refseq
             && self -> num_app_vols [ appREFSEQ ] != 0
             && self -> num_app_vols [ appSRA ] == 0 )
        {
            rc = VResolverForceRemoteRefseq ( self );
        }
    }

    if ( rc == 0 && self -> ticket != NULL && ! have_remote_protected )
        rc = VResolverForceRemoteProtected ( self );

    self -> protocols = eProtocolHttp;

    return rc;
}

#if 1
LIB_EXPORT
rc_t CC VResolverProtocols ( VResolver * self,
    VRemoteProtocols protocols )
{
    if ( self == NULL ) {
        return RC ( rcVFS, rcResolver, rcUpdating, rcSelf, rcNull );
    }

    if ( protocols >= eProtocolLastDefined ) {
        return RC ( rcVFS, rcResolver, rcUpdating, rcParam, rcInvalid );
    }

    self -> protocols = protocols;

    return 0;
}
#endif


/* Make
 *  internal factory function
 */
static
rc_t VResolverMake ( VResolver ** objp, const KDirectory *wd,
    const KRepository *protected, const KConfig *kfg )
{
    rc_t rc;

    VResolver *obj = calloc ( 1, sizeof * obj );
    if ( obj == NULL )
        rc = RC ( rcVFS, rcMgr, rcCreating, rcMemory, rcExhausted );
    else
    {
        VectorInit ( & obj -> roots, 0, 8 );
        VectorInit ( & obj -> local, 0, 8 );
        VectorInit ( & obj -> remote, 0, 8 );
        obj -> wd = wd;

        KRefcountInit ( & obj -> refcount, 1, "VResolver", "make", "resolver" );

        rc = VResolverLoad ( obj, protected, kfg );
        if ( rc == 0 )
        {
            * objp = obj;
            return 0;
        }

        VResolverWhack ( obj );
    }

    return rc;
}

/* Make
 *  ask the VFS manager or repository to make a resolver
 */
LIB_EXPORT
rc_t CC VFSManagerMakeResolver ( const VFSManager * self,
    VResolver ** new_resolver, const KConfig * cfg )
{
    rc_t rc;

    if ( new_resolver == NULL )
        rc = RC ( rcVFS, rcMgr, rcCreating, rcParam, rcNull );
    else
    {
        if ( self == NULL )
            rc = RC ( rcVFS, rcMgr, rcCreating, rcSelf, rcNull );
        else if ( cfg == NULL )
            rc = RC ( rcVFS, rcMgr, rcCreating, rcParam, rcNull );
        else
        {
            KDirectory *wd;
            rc = VFSManagerGetCWD ( self, & wd );
            if ( rc == 0 )
            {
                const KRepositoryMgr *rmgr;
                rc = KConfigMakeRepositoryMgrRead ( cfg, & rmgr );
                if ( rc == 0 )
                {
                    const KRepository *protected = NULL;
                    rc = KRepositoryMgrCurrentProtectedRepository ( rmgr, & protected );
                    if ( rc == 0 || GetRCState ( rc ) == rcNotFound )
                    {
                        rc = VResolverMake ( new_resolver, wd, protected, cfg );
                        KRepositoryRelease ( protected );

                        if ( rc == 0 )
                        {
                            KRepositoryMgrRelease ( rmgr );
                            return 0;
                        }
                    }

                    KRepositoryMgrRelease ( rmgr );
                }

                KDirectoryRelease ( wd );
            }
        }

        *new_resolver = NULL;
    }

    return rc;
}

LIB_EXPORT
rc_t CC KRepositoryMakeResolver ( const KRepository *self,
    VResolver ** new_resolver, const KConfig * cfg )
{
    rc_t rc;

    if ( new_resolver == NULL )
        rc = RC ( rcVFS, rcMgr, rcCreating, rcParam, rcNull );
    else
    {
        if ( self == NULL )
            rc = RC ( rcVFS, rcMgr, rcCreating, rcSelf, rcNull );
        else if ( cfg == NULL )
            rc = RC ( rcVFS, rcMgr, rcCreating, rcParam, rcNull );
        else
        {
            KDirectory *wd;
            rc = KDirectoryNativeDir ( & wd );
            if ( rc == 0 )
            {
                rc = VResolverMake ( new_resolver, wd, self, cfg );
                if ( rc == 0 )
                    return 0;

                KDirectoryRelease ( wd );
            }
        }

        *new_resolver = NULL;
    }

    return rc;
}
